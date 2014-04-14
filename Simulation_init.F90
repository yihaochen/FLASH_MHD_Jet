!!****if* source/Simulation/SimulationMain/WindTunnel/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init(integer(in) :: myPE)
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the wind tunnel with a step problem
!!
!! ARGUMENTS
!!
!!  myPE - the local processor number
!!
!! PARAMETERS
!!
!!  sim_pAmbient    Initial ambient pressure
!!  sim_rhoAmbient  Initial ambient density
!!  sim_windVel     Inflow velocity (parallel to x-axis)
!!  gamma           the Gamma EOS thing
!!  smallp          minimum for pressure
!!  smallx          minimum for abundance
!!
!!***

subroutine Simulation_init(myPE)

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
#include "constants.h"
#include "Flash.h"

  integer :: nozzle = 1
  integer, intent(in) :: myPE

  call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('sim_windVel', sim_windVel)
  call RuntimeParameters_get('sim_pJet', sim(nozzle)%pressure)
  call RuntimeParameters_get('sim_rhoJet', sim(nozzle)%density)
  call RuntimeParameters_get('sim_velJet', sim(nozzle)%velocity)
  call RuntimeParameters_get('sim_nozzleRadius', sim(nozzle)%radius)
  call RuntimeParameters_get('sim_nozzleLength', sim(nozzle)%length)
  call RuntimeParameters_get('sim_nozzlePosX', sim(nozzle)%pos(1))
  call RuntimeParameters_get('sim_nozzlePosY', sim(nozzle)%pos(2))
  call RuntimeParameters_get('sim_nozzlePosZ', sim(nozzle)%pos(3))
  call RuntimeParameters_get('sim_nozzleVecX', sim(nozzle)%vec(1))
  call RuntimeParameters_get('sim_nozzleVecY', sim(nozzle)%vec(2))
  call RuntimeParameters_get('sim_nozzleVecZ', sim(nozzle)%vec(3))
  sim(nozzle)%vec = sim(nozzle)%vec/ sqrt(sum(sim(nozzle)%vec*sim(nozzle)%vec))
  call RuntimeParameters_get('sim_gamma', sim_gamma)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('accx', sim_accx)
  call RuntimeParameters_get('accy', sim_accy)



end subroutine Simulation_init
