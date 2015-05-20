!!****if* source/Simulation/SimulationMain/MHD_Jet/Simulation_initSpecies
!!
!! NAME
!!
!!  Simulation_initSpecies
!!
!!
!! SYNOPSIS
!!  Simulation_initSpecies()
!!
!! DESCRIPTION
!!
!!  This routine will initialize the species and species values needed for
!!  the TwoGamma setup, which advects two fluids with different Gamma values
!!
!!***

subroutine Simulation_initSpecies()

#include "Flash.h"
#include "Multispecies.h"

  use Multispecies_interface
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Simulation_data

  implicit none

  integer :: nozzle=1 

  call RuntimeParameters_get('sim_mu', sim_mu)
  call RuntimeParameters_get('sim_gammaICM', sim_gamma)
  call RuntimeParameters_get('sim_gammaJet', sim(nozzle)%gamma)

  call Multispecies_setProperty(ISM_SPEC, A, sim_mu)
  call Multispecies_setProperty(ISM_SPEC, Z, 1.)
  call Multispecies_setProperty(ISM_SPEC, GAMMA, sim_gamma)

  call Multispecies_setProperty(JET_SPEC, A, sim_mu)
  call Multispecies_setProperty(JET_SPEC, Z, 1.)
  call Multispecies_setProperty(JET_SPEC, GAMMA, sim(nozzle)%gamma)

end subroutine Simulation_initSpecies
