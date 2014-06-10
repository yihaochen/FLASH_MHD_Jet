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

  function taper(nozzle, r, var_in, var_out)

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: r, var_in, var_out
    real :: r1, r2, rjet
    real :: taper
    rjet = sim(nozzle)%radius
    r1 = sim(nozzle)%bfeather_inner
    r2 = rjet - sim(nozzle)%bfeather_outer

    if (r.ge.0.0 .and. r.lt.r1) then
      taper = (-2*r**3/r1**3 + 3*r**2/r1**2)*var_in
    else if (r.ge.r1 .and. r.lt.r2) then
      taper = 1.0*var_in
    else if (r.ge.r2 .and. r.lt.rjet) then
      taper = (-2*(rjet-r)**3/(rjet-r2)**3 &
      + 3*(rjet-r)**2/(rjet-r2)**2)*(var_in-var_out) + var_out
    else
      taper = var_out
    end if

  end function taper

end module Simulation_data


