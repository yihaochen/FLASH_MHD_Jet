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
      real :: bfeather_inner, bfeather_outer, zfeather
      real :: bphi, bz
  end TYPE nozzle_struct

  !! *** Runtime Parameters *** !!

  type(nozzle_struct), save, dimension(NOZZLES) :: sim

  real, save :: sim_pAmbient, sim_rhoAmbient, sim_windVel, sim_bzAmbient
  real, save :: sim_gamma, sim_smallP, sim_smallX

contains

  function cross(a,b)

    real, dimension(3), INTENT(in) :: a,b
    real, dimension(3) :: cross
    cross(:)=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]

  end function cross

  function taperR(nozzle, r, var_in, var_out)

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: r, var_in, var_out
    real :: r1, r2, rout
    real :: taperR
    r1 = sim(nozzle)%bfeather_inner
    r2 = sim(nozzle)%radius
    rout = sim(nozzle)%radius + sim(nozzle)%bfeather_outer

    if (r.ge.0.0 .and. r.lt.r1) then
      taperR = (-2*r**3/r1**3 + 3*r**2/r1**2)*var_in
    else if (r.ge.r1 .and. r.lt.r2) then
      taperR = 1.0*var_in
    else if (r.ge.r2 .and. r.lt.rout) then
      taperR = (-2*(rout-r)**3/(rout-r2)**3 &
      + 3*(rout-r)**2/(rout-r2)**2)*(var_in-var_out) + var_out
    else
      taperR = var_out
    end if

  end function taperR

  function taper(nozzle, r, z, var_cen, var_in, var_out)
  ! Taper function for both R and z direction
  ! For radial direction, variable could have 3 values: center, inside, and outside
  !
  ! 0          bfeather_inner       radius         radius+bfeather_outer
  ! 0               r1                r2                   rout
  ! |   center       |       inside    |    outside          |
  ! |   var_cen      |       var_in    |    var_out          |

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: r, z, var_cen, var_in, var_out
    real :: r1, r2, rout, zout, zjet, var_zcen, var_zin
    real :: taper
    r1 = sim(nozzle)%bfeather_inner
    r2 = sim(nozzle)%radius 
    rout = sim(nozzle)%radius + sim(nozzle)%bfeather_outer
    zjet = sim(nozzle)%length
    zout = zjet + sim(nozzle)%zfeather

    ! z part
    if (abs(z).ge.0.0 .and. abs(z).lt.zjet) then
      var_zin = 1.0*var_in
      var_zcen = 1.0*var_cen
    else if (abs(z).ge.zjet .and. abs(z).lt.zout) then
      var_zin = (-2*(zout-abs(z))**3/(zout-zjet)**3 &
      + 3*(zout-abs(z))**2/(zout-zjet)**2)*(var_in-var_out) + var_out
      var_zcen = (-2*(zout-abs(z))**3/(zout-zjet)**3 &
      + 3*(zout-abs(z))**2/(zout-zjet)**2)*(var_cen-var_out) + var_out
    else
      var_zin = var_out
      var_zcen = var_out
    end if

    ! radial part
    if (r.ge.0.0 .and. r.lt.r1) then
      taper = (-2*r**3/r1**3 + 3*r**2/r1**2)*(var_zin-var_zcen) + var_zcen
    else if (r.ge.r1 .and. r.lt.r2) then
      taper = 1.0*var_zin
    else if (r.ge.r2 .and. r.lt.rout) then
      taper = (-2*(rout-r)**3/(rout-r2)**3 &
      + 3*(rout-r)**2/(rout-r2)**2)*(var_zin-var_out) + var_out
    else
      taper = var_out
    endif


  end function taper

end module Simulation_data


