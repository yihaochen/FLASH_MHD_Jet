!!****if* source/Simulation/SimulationMain/WindTunnel/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for the wind tunnel problem with a step
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
!!
!!   
!!
!!***

module Simulation_data

  implicit none
#include "Simulation.h"

  TYPE nozzle_struct
      real, dimension(3) :: pos, jetvec
      real :: radius, length
      real :: pressure, density, velocity, gamma
      real :: bfeather_inner, bfeather_outer
      real :: bphi
  end TYPE nozzle_struct

  !! *** Runtime Parameters *** !!

  type(nozzle_struct), save, dimension(NOZZLES) :: sim

  real, save :: sim_pAmbient, sim_rhoAmbient, sim_windVel
  !real, save :: sim_pJet, sim_rhoJet, sim_velJet
  real, save :: sim_gamma, sim_smallP, sim_smallX
  
contains
  
  function cross(a,b)
    
    real, dimension(3), INTENT(in) :: a,b
    real, dimension(3) :: cross
    cross(:)=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
    
  end function cross
  

end module Simulation_data


