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

    implicit none

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: time, dt
    real :: p, g, v, r, L, bf

    ! --------------------------------------------------------------------------
    ! Update the hydro variables of the jet nozzle according to the wind-driven bubble solution.
    g = sim(nozzle)%gamma
    v = sim(nozzle)%velocity
    r = sim(nozzle)%radius
    bf= sim(nozzle)%rFeatherOut
    L = sim(nozzle)%power

    if (sim(nozzle)%t0 < 0.0) then
       sim(nozzle)%t0 = ((r+0.5*bf)**2*PI*2*v)**1.25*(sim_rhoAmbient*g/(g-1)/L)**0.75*0.227082*2.0
    endif
    ! Calculate the jet pressure
    sim(nozzle)%pressure = (max(time,sim(nozzle)%t0))**(-0.8)&
                           *0.305454*sim_rhoAmbient**0.6*((g-1)/g*L)**0.4
    sim(nozzle)%density = 2.0/v/v*(L/((r+0.5*bf)**2*PI*2*v) - g/(g-1)*sim(nozzle)%pressure)
    !sim(nozzle)%density = max(gr_smallrho, sim(nozzle)%density)
    sim(nozzle)%mach = v/sqrt(g*sim(nozzle)%pressure/sim(nozzle)%density)

    sim(nozzle)%bzOld = sim(nozzle)%bz
    sim(nozzle)%bz = sqrt(2.0*sim(nozzle)%pressure/sim(nozzle)%beta/(1.0+sim(nozzle)%helicity**2))
    sim(nozzle)%bphi = sim(nozzle)%bz*sim(nozzle)%helicity
    !sim(nozzle)%deltaP = ( (max(time,t0))**(-0.8)-(max(time-dt,t0))**(-0.8) )&
    !                       *0.305454*sim_rhoAmbient**0.6*((g-1)/g*L)**0.4
    !sim(nozzle)%deltaRho = -2.0/v/v*(g/(g-1)*sim(nozzle)%deltaP)


    ! --------------------------------------------------------------------------
    ! Update the jet nozzle position and direction according to the velocity and
    ! angular velocity.
    sim(nozzle)%posOld = sim(nozzle)%pos
    sim(nozzle)%jetvecOld = sim(nozzle)%jetvec

    sim(nozzle)%pos = sim(nozzle)%pos + sim(nozzle)%linVel*dt
    sim(nozzle)%jetvec = sim(nozzle)%jetvec + cross(sim(nozzle)%angVel, sim(nozzle)%jetvec)*dt

  end subroutine sim_jetNozzleUpdate

end module Simulation_jetNozzleUpdate
