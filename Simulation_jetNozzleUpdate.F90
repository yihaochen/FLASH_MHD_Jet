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
    use Driver_data, ONLY: dr_dtMin

    implicit none

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: time, dt
    ! Calculate the jet pressure
    real :: p, g, v, r, L
    !real, dimension(3) :: angVel = (/ 0.0, 0.0, 2.0*PI/1.0E13 /)
    ! Angular velocity
    sim(nozzle)%angVel = (/ 0.0, 0.0, 0.0 /)
    ! Linear velocity
    sim(nozzle)%linVel = (/ 0.0, 6.0E8, 0.0 /)

    ! --------------------------------------------------------------------------
    ! Update the hydro variables of the jet nozzle according to the wind-driven bubble solution.
    g = sim(nozzle)%gamma
    v = sim(nozzle)%velocity
    r = sim(nozzle)%radius
    L = sim(nozzle)%power

    if (t0 < 0.0) then
       t0 = (r*r*PI*2*v)**1.25*(sim_rhoAmbient*g/(g-1)/L)**0.75*0.227082*2.0
    endif
    sim(nozzle)%pressure = (max(time,t0))**(-0.8)*0.305454*sim_rhoAmbient**0.6*((g-1)/g*L)**0.4
    sim(nozzle)%density = 2.0/v/v*(L/(r*r*PI*2*v) - g/(g-1)*sim(nozzle)%pressure)
    !sim(nozzle)%density = max(gr_smallrho, sim(nozzle)%density)
    sim(nozzle)%mach = v/sqrt(g*sim(nozzle)%pressure/sim(nozzle)%density)

    sim(nozzle)%deltaP = ( (max(time,t0))**(-0.8)-(max(time-dt,t0))**(-0.8) )&
                           *0.305454*sim_rhoAmbient**0.6*((g-1)/g*L)**0.4
    sim(nozzle)%deltaRho = -2.0/v/v*(g/(g-1)*sim(nozzle)%deltaP)


    ! --------------------------------------------------------------------------
    ! Update the jet nozzle position and direction according to the velocity and
    ! angular velocity.
    sim(nozzle)%pos_old = sim(nozzle)%pos
    sim(nozzle)%jetvec_old = sim(nozzle)%jetvec

    sim(nozzle)%pos = sim(nozzle)%pos + sim(nozzle)%linVel*dt
    sim(nozzle)%jetvec = sim(nozzle)%jetvec + cross(sim(nozzle)%angVel, sim(nozzle)%jetvec)*dt

  end subroutine sim_jetNozzleUpdate

end module Simulation_jetNozzleUpdate
