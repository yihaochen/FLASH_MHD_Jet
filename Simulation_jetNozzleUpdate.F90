!!****if* source/Simulation/SimulationMain/magnetoHD/MHD_Jet/Simulation_jetNozzleUpdate
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
    use Driver_data, ONLY : dr_nStep, dr_restart
    use Hydro_data, ONLY : hy_bref

    implicit none

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: time, dt
    real :: p, g, v, R, L, bf, M, dt_start, dt_end, h, x

    integer :: funit = 99, isFirst = 1
    integer :: ioStat

    character (len=13), save :: nozzleVecFName
    nozzleVecFName = 'nozzleVec.dat'

    ! --------------------------------------------------------------------------
    ! Update the hydro variables of the jet nozzle
    if (time.ge.sim(nozzle)%tOn .and. time.lt.sim(nozzle)%tOn+sim(nozzle)%duration) then
       sim(nozzle)%on = .true.
    else
       sim(nozzle)%on = .false.
    endif

    dt_start = sim(nozzle)%duration/20.0
    dt_end = sim(nozzle)%duration/100.0
    g = sim(nozzle)%gamma
    ! Target velocity of the jet
    v = sim(nozzle)%velJet
    R = sim(nozzle)%radius
    bf= sim(nozzle)%rFeatherOut
    L = sim(nozzle)%power

    h = sim(nozzle)%helicity
    ! Ratio of (enthalpy + toroidal magnetic energy density) to pressure
    x = g/(g-1.0)+h**2/(1.0+h**2)/sim(nozzle)%beta


    if (sim(nozzle)%on) then
       sim(nozzle)%velocity = sim(nozzle)%velJet
       ! Turn on the jet by increasing the mach number from initMach to mach.
       if (time .lt. sim(nozzle)%tOn+dt_start) then
          !M = sim(nozzle)%mach
          !M = sim(nozzle)%initMach + (sim(nozzle)%mach-sim(nozzle)%initMach)&
          !    *cos(PI*( max(-0.5, min(0.0, 0.5*(time-sim(nozzle)%tOn-t1)/t1))))
          M = max(sim(nozzle)%initMach, sim(nozzle)%mach*((time-sim(nozzle)%tOn)/dt_start)**0.7)
       else
          M = sim(nozzle)%mach
       endif

       ! The density and pressure are calculated using the target velocity,
       ! so that they will remain constant during jet turning off time.
       sim(nozzle)%density = 0.5*L/PI/v**3/( R*R*(0.5+x/M**2/g) + R*bf*(0.3125+x/M**2/g) &
                             + bf*bf*(0.06056+0.29736*x/M**2/g) )
       sim(nozzle)%pressure = v*v*sim(nozzle)%density/M**2/g

       ! Turn off the jet by decreasing the velocity
       if (time .gt. sim(nozzle)%tOn+sim(nozzle)%duration-dt_end) then
           sim(nozzle)%velocity = sim(nozzle)%velJet&
           *cos(PI*( max(0.0, min(0.5, 0.5*(time-sim(nozzle)%tOn-sim(nozzle)%duration+dt_end)/dt_end))))
           if (sim_meshMe  == MASTER_PE .and. time .lt. sim(nozzle)%tOn+sim(nozzle)%duration) then
               write(*,*) '[Jet turning off] current v = ', sim(nozzle)%velocity
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
       sim(nozzle)%bz = sqrt(2.0*sim(nozzle)%pressure/sim(nozzle)%beta/(1.0+sim(nozzle)%helicity**2))*hy_bref
       sim(nozzle)%bphi = sim(nozzle)%bz*sim(nozzle)%helicity

       ! --------------------------------------------------------------------------
       ! Update the jet nozzle position and direction according to the velocity and
       ! angular velocity.
       !sim(nozzle)%jetvecOld = sim(nozzle)%jetvec

       if (dt.gt.0.0 .and. sim(nozzle)%precangle.gt.0.0) then
          if (sim_useTableJiggle) then
             call Timers_start('Simulation_jiggleRead')
             call Simulation_jiggleRead(nozzle, time, dt)
             call Timers_stop('Simulation_jiggleRead')
          else
             call Timers_start('Simulation_jiggle')
             call Simulation_jiggle(nozzle, time, dt)
             call Timers_stop('Simulation_jiggle')
          endif
          sim(nozzle)%posOld = sim(nozzle)%pos
          sim(nozzle)%pos = sim(nozzle)%pos + sim(nozzle)%linVel*dt
          !sim(nozzle)%jetvec = sim(nozzle)%jetvec + cross(sim(nozzle)%angVel, sim(nozzle)%jetvec)*dt

          ! Write the jet nozzle vectors to file.
          if (sim_meshMe == MASTER_PE) then

             ! create the file from scratch if it is a not a restart simulation,
             ! otherwise append to the end of the file

             !No mater what, we are opening the file. Check to see if already there
             ioStat = 0
             open(funit, file=trim(nozzleVecFName), position='APPEND', status='OLD', iostat=ioStat)
             if(ioStat .NE. 0) then
                !print *, 'FILE FOUND'
                open(funit, file=trim(nozzleVecFName), position='APPEND')
             endif

             ! Header
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
10              format (2X, 1(a10, :, 1X), 7(a25, :, 1X), 1(a10, :, 1X))
                isFirst = 0

             ! Restart mark
             else if(isFirst .EQ. 1) then
                write (funit, 11)
11              format('# simulation restarted')
                isFirst = 0
             endif

             ! Actual write
             write (funit, 12) dr_nstep, time, sim(nozzle)%jetvec, sim(nozzle)%angVel, sim(nozzle)%randSeed
12           format (1X, 1(I10, :, 1X), 7(es25.18, :, 1X), 1(I10, :, 1X))

             close (funit)          ! Close the file.
          endif !Write nozzle vector to file
       endif !Nozzle jiggle is on
   endif !Nozzle is on

  end subroutine sim_jetNozzleUpdate

end module Simulation_jetNozzleUpdate
