!!****if* source/Simulation/SimulationMain/MHD_Jet/Simulation_jetNozzleUpdate
!!
!! NAME
!!  Simulation_jetNozzleUpdate
!!
!! SYNOPSIS
!!
!!  use Simulation_jetNozzleUpdate
!!
!! DESCRIPTION
!!
!!  Handle the movement and the rotation of the nozzle.
!!
!! ARGUMENTS
!!
!!
!!
!!
!!
!!***

module Simulation_jetNozzleUpdate

contains

  subroutine sim_jetNozzleUpdate(nozzle, time, dt)

#include "constants.h"
#include "Simulation.h"

    use Simulation_data
    use Simulation_interface, ONLY: Simulation_jiggle
    use Timers_interface, ONLY : Timers_start, Timers_stop
    use Driver_data, ONLY : dr_globalMe, dr_nStep, dr_restart

    implicit none

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: time, dt
    real, dimension(3) :: nutationVec
    real :: p, g, v, r, L, bf, M, t1

    integer :: funit = 99, isFirst = 1
    integer :: ioStat

    character (len=13), save :: nozzleVecFName
    nozzleVecFName = 'nozzleVec.dat'

    ! --------------------------------------------------------------------------
    ! Update the hydro variables of the jet nozzle (according to the wind-driven bubble solution.)

    !if (sim(nozzle)%density < 0.0) then
       !sim(nozzle)%t0 = (r**2*PI*2*v)**1.25*(sim_rhoCore*g/(g-1)/L)**0.75*0.227082


    t1 = sim(nozzle)%duration/100.0
    ! Turn on the jet by increasing the mach number from initMach to mach.
    if (time .lt. sim(nozzle)%tOn+t1) then
       g = sim(nozzle)%gamma
       v = sim(nozzle)%velocity
       r = sim(nozzle)%radius
       bf= sim(nozzle)%rFeatherOut
       L = sim(nozzle)%power

       !M = sim(nozzle)%mach
       M = sim(nozzle)%initMach + (sim(nozzle)%mach-sim(nozzle)%initMach)&
           *cos(PI*( max(-0.5, min(0.0, 0.5*(time-sim(nozzle)%tOn-t1)/t1))))

       sim(nozzle)%density = 0.5*L/PI/v**3/( 0.5*r*r*(1.+1./M**2/(g-1.)) + r*bf*(0.3125+1./M**2/(g-1.)) &
                             + bf*bf*(0.06056+0.29736/M**2/(g-1.)) )
       sim(nozzle)%pressure = v*v*sim(nozzle)%density/M**2/g
    endif

    ! Turn off the jet by decreasing the velocity 
    if (time .gt. sim(nozzle)%tOn+sim(nozzle)%duration-t1) then
        sim(nozzle)%velocity = sim(nozzle)%velJet&
        *cos(PI*( max(0.0, min(0.5, 0.5*(time-sim(nozzle)%tOn-sim(nozzle)%duration+t1)/t1))))
        if (dr_globalMe  == MASTER_PE .and. time .lt. sim(nozzle)%tOn+sim(nozzle)%duration) then
            write(*,*) 'v = ', sim(nozzle)%velocity
        endif


    endif
    !endif
    ! Calculate the jet pressure


    !sim(nozzle)%pressure = (max(time,sim(nozzle)%t0))**(-0.8)&
    !                       *0.305454*sim_rhoCore**0.6*((g-1)/g*L)**0.4
    !sim(nozzle)%density = 2.0/v/v*(L/(r**2*PI*2*v) - g/(g-1)*sim(nozzle)%pressure)
    !sim(nozzle)%density = max(gr_smallrho, sim(nozzle)%density)
    !sim(nozzle)%mach = v/sqrt(g*sim(nozzle)%pressure/sim(nozzle)%density)

    !sim(nozzle)%deltaP = ( (max(time,t0))**(-0.8)-(max(time-dt,t0))**(-0.8) )&
    !                       *0.305454*sim_rhoCore**0.6*((g-1)/g*L)**0.4
    !sim(nozzle)%deltaRho = -2.0/v/v*(g/(g-1)*sim(nozzle)%deltaP)

    sim(nozzle)%bzOld = sim(nozzle)%bz
    sim(nozzle)%bz = sqrt(2.0*sim(nozzle)%pressure/sim(nozzle)%beta/(1.0+sim(nozzle)%helicity**2))
    sim(nozzle)%bphi = sim(nozzle)%bz*sim(nozzle)%helicity


    ! --------------------------------------------------------------------------
    ! Update the jet nozzle position and direction according to the velocity and
    ! angular velocity.
    !sim(nozzle)%jetvecOld = sim(nozzle)%jetvec

    if (dt.gt.0.0) then
       call Timers_start('Simulation_jiggle')
       call Simulation_jiggle(nozzle,  time, dt)
       call Timers_stop('Simulation_jiggle')
       sim(nozzle)%posOld = sim(nozzle)%pos
       sim(nozzle)%pos = sim(nozzle)%pos + sim(nozzle)%linVel*dt
       !sim(nozzle)%jetvec = sim(nozzle)%jetvec + cross(sim(nozzle)%angVel, sim(nozzle)%jetvec)*dt

       ! Write the jet nozzle vectors to file.
       if (dr_globalMe  == MASTER_PE) then

          ! create the file from scratch if it is a not a restart simulation,
          ! otherwise append to the end of the file

          !No mater what, we are opening the file. Check to see if already there
          ioStat = 0
          open(funit, file=trim(nozzleVecFName), position='APPEND', status='OLD', iostat=ioStat)
          if(ioStat .NE. 0) then
             !print *, 'FILE FOUND'
             open(funit, file=trim(nozzleVecFName), position='APPEND')
          endif

          if (isFirst .EQ. 1 .AND. (.NOT. dr_restart .or. ioStat .NE. 0)) then
             write (funit, 10)   &
                  '#step      ', &
                  'time                     ', &
                  'nozzleVecX               ', &
                  'nozzleVecY               ', &
                  'nozzleVecZ               ', &
                  'nozzleAngVelX            ', &
                  'nozzleAngVelY            ', &
                  'nozzleAngVelZ            ', &
                  'randSeed  '
10           format (2X, 1(a10, :, 1X), 7(a25, :, 1X), 1(a10, :, 1X))

          else if(isFirst .EQ. 1) then
             write (funit, 11)
11           format('# simulation restarted')
          endif

          ! Write the global sums to the file.
          write (funit, 12) dr_nstep, time, sim(nozzle)%jetvec, sim(nozzle)%angVel, sim(nozzle)%randSeed
12        format (1X, 1(I10, :, 1X), 7(es25.18, :, 1X), 1(I10, :, 1X))

          close (funit)          ! Close the file.
          isFirst = 0
       endif
    endif


  end subroutine sim_jetNozzleUpdate

end module Simulation_jetNozzleUpdate
