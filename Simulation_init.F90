!!****if* source/Simulation/SimulationMain/MHD_Jet/Simulation_init
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
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Simulation_jetNozzleUpdate, ONLY : sim_jetNozzleUpdate
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_data, ONLY : dr_globalMe, dr_simTime, dr_dt, dr_restart
  use Grid_data, ONLY : gr_minCellSize
  use IO_interface, ONLY :  IO_getScalar
  !use Grid_data, ONLY : gr_smallrho

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Simulation.h"

  integer :: nozzle = 1


  call RuntimeParameters_get('smlrho', sim_smlrho)
  call RuntimeParameters_get('smallp', sim_smallp)
  call RuntimeParameters_get('smalle', sim_smalle)
  call RuntimeParameters_get('smallx', sim_smallx)

  call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('sim_windVel', sim_windVel)
  call RuntimeParameters_get('sim_gammaAmbient', sim_gamma)
  call RuntimeParameters_get('sim_bzAmbient', sim_bzAmbient)
  call RuntimeParameters_get('sim_densityProfile', sim_densityProfile)
  if (sim_densityProfile == "betacore") then
     call RuntimeParameters_get('sim_densityBeta', sim_densityBeta)
     call RuntimeParameters_get('sim_densityCoreR', sim_densityCoreR)
  endif
  call RuntimeParameters_get('sim_powerJet', sim(nozzle)%power)
  !call RuntimeParameters_get('sim_rhoJet', sim(nozzle)%density)
  call RuntimeParameters_get('sim_velJet', sim(nozzle)%velocity)
  call RuntimeParameters_get('sim_machJet', sim(nozzle)%mach)
  call RuntimeParameters_get('sim_outflowRatio', sim(nozzle)%outflowR)
  call RuntimeParameters_get('sim_gammaJet', sim(nozzle)%gamma)
  call RuntimeParameters_get('sim_betaJet', sim(nozzle)%beta)
  call RuntimeParameters_get('sim_helicityJet', sim(nozzle)%helicity)
  call RuntimeParameters_get('sim_timeMHDon', sim(nozzle)%timeMHDon)

  call RuntimeParameters_get('nozzleRadius', sim(nozzle)%radius)
  call RuntimeParameters_get('nozzleHalfL', sim(nozzle)%length)
  call RuntimeParameters_get('zTorInj', sim(nozzle)%zTorInj)
  call RuntimeParameters_get('rFeatherIn', sim(nozzle)%rFeatherIn)
  call RuntimeParameters_get('rFeatherOut', sim(nozzle)%rFeatherOut)
  sim(nozzle)%rFeatherMix = gr_minCellSize*NGUARD
  call RuntimeParameters_get('zFeather', sim(nozzle)%zFeather)
  sim(nozzle)%zFeatherMix = gr_minCellSize
  call RuntimeParameters_get('initGeometry', sim(nozzle)%initGeometry)
  call RuntimeParameters_get('derefine_z1', sim(nozzle)%derefine_z1)
  call RuntimeParameters_get('derefine_z2', sim(nozzle)%derefine_z2)
  call RuntimeParameters_get('refine_jetR', sim(nozzle)%refine_jetR)
  call RuntimeParameters_get('derefine_jetR', sim(nozzle)%derefine_jetR)
  call RuntimeParameters_get('lrefine_0', sim(nozzle)%lrefine_0)
  call RuntimeParameters_get('sim_ptInitNum', sim_ptInitNum)
  call RuntimeParameters_get('sim_ptAddPeriod', sim_ptAddPeriod)

  if (dr_restart) then
     call IO_getScalar('nozzlePosX', sim(nozzle)%pos(1))
     call IO_getScalar('nozzlePosY', sim(nozzle)%pos(2))
     call IO_getScalar('nozzlePosZ', sim(nozzle)%pos(3))
     call IO_getScalar('nozzleVecX', sim(nozzle)%jetvec(1))
     call IO_getScalar('nozzleVecY', sim(nozzle)%jetvec(2))
     call IO_getScalar('nozzleVecZ', sim(nozzle)%jetvec(3))
     call IO_getScalar('nozzleAngVelX', sim(nozzle)%angVel(1))
     call IO_getScalar('nozzleAngVelY', sim(nozzle)%angVel(2))
     call IO_getScalar('nozzleAngVelZ', sim(nozzle)%angVel(3))
     call IO_getScalar('nozzleLinVelX', sim(nozzle)%linVel(1))
     call IO_getScalar('nozzleLinVelY', sim(nozzle)%linVel(2))
     call IO_getScalar('nozzleLinVelZ', sim(nozzle)%linVel(3))
     call IO_getScalar('nozzlePressure', sim(nozzle)%pressure)
     call IO_getScalar('nozzleDensity', sim(nozzle)%density)
     call IO_getScalar('nozzleBz', sim(nozzle)%bz)
     call IO_getScalar('nozzleBphi', sim(nozzle)%bphi)
     !call IO_getScalar('nozzlet0', sim(nozzle)%t0)

  else

     call RuntimeParameters_get('nozzlePosX', sim(nozzle)%pos(1))
     call RuntimeParameters_get('nozzlePosY', sim(nozzle)%pos(2))
     call RuntimeParameters_get('nozzlePosZ', sim(nozzle)%pos(3))
     call RuntimeParameters_get('nozzleVecX', sim(nozzle)%jetvec(1))
     call RuntimeParameters_get('nozzleVecY', sim(nozzle)%jetvec(2))
     call RuntimeParameters_get('nozzleVecZ', sim(nozzle)%jetvec(3))
     sim(nozzle)%jetvec = sim(nozzle)%jetvec/ sqrt(sum(sim(nozzle)%jetvec*sim(nozzle)%jetvec))
     call RuntimeParameters_get('nozzleAngVelX', sim(nozzle)%angVel(1))
     call RuntimeParameters_get('nozzleAngVelY', sim(nozzle)%angVel(2))
     call RuntimeParameters_get('nozzleAngVelZ', sim(nozzle)%angVel(3))
     call RuntimeParameters_get('nozzleLinVelX', sim(nozzle)%linVel(1))
     call RuntimeParameters_get('nozzleLinVelY', sim(nozzle)%linVel(2))
     call RuntimeParameters_get('nozzleLinVelZ', sim(nozzle)%linVel(3))
     sim(nozzle)%density = -1.0
     call sim_jetNozzleUpdate(nozzle, dr_simTime, dr_dt)

     if (sim(nozzle)%zTorInj < sim(nozzle)%length+1.5*sim(nozzle)%zFeather .and.&
        dr_globalMe==MASTER_PE) then
        print*, '!!!!!!!!'
        print*, 'Warning! zTorInj is too small that it overlaps with the nozzle.'
        print*, 'Toroidal field will be smaller than it should be.'
        print*, '!!!!!!!!'
     endif
  endif

  if (dr_globalMe==MASTER_PE) then
     !write(*,'(a, es11.3)') 't0:', sim(nozzle)%t0
     write(*,'(a, 2es11.3, f7.2)') '(p, rho, M)=', &
     sim(nozzle)%pressure, sim(nozzle)%density, sim(nozzle)%mach
     write(*,'(a, 2es11.3)') '(bz, bphi)=', sim(nozzle)%bz, sim(nozzle)%bphi
     write(*,'(a, es11.3)') 'rFeatherOut:' , sim(nozzle)%rFeatherOut
     write(*,'(a, es11.3)') 'rFeatherMix:' , sim(nozzle)%rFeatherMix
     write(*,'(a, es11.3)') 'zFeatherMix:' , sim(nozzle)%zFeatherMix
  endif


end subroutine Simulation_init
