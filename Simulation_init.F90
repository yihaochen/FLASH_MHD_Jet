!!****if* source/Simulation/SimulationMain/WindTunnel/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the wind tunnel with a step problem
!!
!! ARGUMENTS
!!
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

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_data, ONLY : dr_globalMe
  !use Grid_data, ONLY : gr_smallrho

  implicit none
#include "constants.h"
#include "Flash.h"

  integer :: nozzle = 1

  call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('sim_windVel', sim_windVel)
  call RuntimeParameters_get('sim_gammaAmbient', sim_gamma)
  call RuntimeParameters_get('sim_bzAmbient', sim_bzAmbient)

  call RuntimeParameters_get('sim_powerJet', sim(nozzle)%power)
  !call RuntimeParameters_get('sim_rhoJet', sim(nozzle)%density)
  call RuntimeParameters_get('sim_velJet', sim(nozzle)%velocity)
  call RuntimeParameters_get('sim_gammaJet', sim(nozzle)%gamma)
  call RuntimeParameters_get('sim_bphiJet', sim(nozzle)%bphi)
  call RuntimeParameters_get('sim_bzJet', sim(nozzle)%bz)
  call RuntimeParameters_get('sim_timeMHDon', sim(nozzle)%timeMHDon)

  call RuntimeParameters_get('nozzleRadius', sim(nozzle)%radius)
  call RuntimeParameters_get('nozzleHalfL', sim(nozzle)%length)
  call RuntimeParameters_get('nozzlePosX', sim(nozzle)%pos(1))
  call RuntimeParameters_get('nozzlePosY', sim(nozzle)%pos(2))
  call RuntimeParameters_get('nozzlePosZ', sim(nozzle)%pos(3))
  call RuntimeParameters_get('nozzleVecX', sim(nozzle)%jetvec(1))
  call RuntimeParameters_get('nozzleVecY', sim(nozzle)%jetvec(2))
  call RuntimeParameters_get('nozzleVecZ', sim(nozzle)%jetvec(3))
  sim(nozzle)%jetvec = sim(nozzle)%jetvec/ sqrt(sum(sim(nozzle)%jetvec*sim(nozzle)%jetvec))
  call RuntimeParameters_get('bPosZ', sim(nozzle)%bPosZ)
  call RuntimeParameters_get('bfeather_inner', sim(nozzle)%bfeather_inner)
  call RuntimeParameters_get('bfeather_outer', sim(nozzle)%bfeather_outer)
  call RuntimeParameters_get('zfeather', sim(nozzle)%zfeather)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smallp', sim_smallP)

  if (sim(nozzle)%bPosZ < sim(nozzle)%length+1.5*sim(nozzle)%zfeather .and.&
     dr_globalMe==MASTER_PE) then
     print*, '!!!!!!!!'
     print*, 'Warning! bPosZ is too small that it overlaps with the nozzle.'
     print*, 'Toroidal field will be smaller than sim_bphiJet'
     print*, '!!!!!!!!'
  endif
  call calc_jet(nozzle, 0.0)
  if (dr_globalMe==MASTER_PE) then
     write(*,'(a, es11.3)') 't0:', t0
     write(*,'(a, 2es11.3)') '(p, rho)=', sim(nozzle)%pressure, sim(nozzle)%density
  endif


end subroutine Simulation_init
