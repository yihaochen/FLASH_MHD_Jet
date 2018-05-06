!!****if* source/Simulation/SimulationMain/MHD_Jet/Simulation_data
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
!!
!!
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"
#include "Simulation.h"

  TYPE nozzle_struct
      logical :: on
      real, dimension(3) :: pos, precvec, jetvec, coneVec 
      real, dimension(3) :: posOld, jetvecOld
      real, dimension(3) :: angVel, linVel
      integer, dimension(1) :: randSeed
      real :: precession, nutation, precangle
      real :: radius, length, zTorInj, outflowR
      !real :: derefine_z1, derefine_z2
      real :: refine_jetR, derefine_jetR
      integer :: lrefine_0
      real :: power, pressure, density, velocity, velJet, gamma, mach, initMach
      real :: rFeatherIn, rFeatherOut, rFeatherMix, zFeather, zFeatherMix
      real :: beta, helicity
      real :: bphi, bz, tOn, duration, timeMHDon
      real :: bzOld
      !real :: t0 = -1.0
      !TODO: initGeometry -> sim_initGeometry
      character(len=MAX_STRING_LENGTH) :: initGeometry

  end TYPE nozzle_struct

  !! *** Runtime Parameters *** !!

  type(nozzle_struct), save, dimension(NOZZLES) :: sim

  real,save :: sim_Tcore, sim_Tout, sim_rhoCore, sim_windVel, sim_bzAmbient
  real,save :: sim_densityBeta, sim_rCore, sim_rCoreT, sim_rCut, sim_mu
  real,save :: sim_gamma, sim_smallp, sim_smallx, sim_smlrho, sim_smalle
  real,save :: sim_rhoFloor
  real,save,allocatable,dimension(:) :: sim_xcoord
  real,save,allocatable,dimension(:) :: sim_ycoord
  real,save,allocatable,dimension(:) :: sim_zcoord
  real,save,allocatable,dimension(:) :: sim_xcoordf
  real,save,allocatable,dimension(:) :: sim_ycoordf
  real,save,allocatable,dimension(:) :: sim_zcoordf
  integer,save :: sim_ptInitNum
  real,save :: sim_ptAddPeriod, sim_ptAddArea, sim_ptMaxRadius
  character(len=MAX_STRING_LENGTH),save :: sim_densityProfile
  logical,save :: sim_useTableJiggle
  integer,save :: sim_nozfileunit=98
  character(len=MAX_STRING_LENGTH),save :: sim_nozVecInput
  integer,save :: sim_meshMe

contains

  function cross(a,b)

    real, dimension(3), INTENT(in) :: a,b
    real, dimension(3) :: cross
    cross(:)=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]

  end function cross

  function taperR(nozzle, r, var_in, var_out, zero_center)

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: r, var_in, var_out
    logical, INTENT(in), optional :: zero_center
    logical :: z_c
    real :: r1, r2, rout, rmix
    real :: taperR

    z_c = .true.
    if (present(zero_center)) z_c = zero_center

    r1 = sim(nozzle)%rFeatherIn
    r2 = sim(nozzle)%radius
    rout = sim(nozzle)%radius + sim(nozzle)%rFeatherOut
    rmix = sim(nozzle)%radius + sim(nozzle)%rFeatherOut + sim(nozzle)%rFeatherMix

    if (r.ge.0.0 .and. r.lt.r1) then
      if (z_c) then
        taperR = (-2*r**3/r1**3 + 3*r**2/r1**2)*var_in
      else
        taperR = 1.0*var_in
      end if
    else if (r.ge.r1 .and. r.lt.rout) then
      taperR = 1.0*var_in
    else if (r.ge.rout .and. r.lt.rmix) then
      taperR = (-2*(rmix-r)**3/(rmix-rout)**3 &
      + 3*(rmix-r)**2/(rmix-rout)**2)*(var_in-var_out) + var_out
    else
      taperR = var_out
    end if

  end function taperR

  function taperL(nozzle, z, var_in, var_out)

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: z, var_in, var_out
    real :: zmix, zjet
    real :: taperL
    zjet = sim(nozzle)%length
    zmix = zjet + sim(nozzle)%zFeatherMix

    ! z part
    if (abs(z).ge.0.0 .and. abs(z).lt.zjet) then
      taperL = 1.0*var_in
    else if (abs(z).ge.zjet .and. abs(z).lt.zmix) then
      taperL = (-2*(zmix-abs(z))**3/(zmix-zjet)**3 &
      + 3*(zmix-abs(z))**2/(zmix-zjet)**2)*(var_in-var_out) + var_out
    else
      taperL = var_out
    end if

  end function taperL
  
  function taper(nozzle, r, z, var_cen, var_in, var_out)
  ! Taper function for both R and z direction
  ! For radial direction, variable could have 3 values: center, inside, and outside
  !
  ! 0          rFeatherIn           radius         radius+rFeatherOut    ..+FeatherMix
  ! 0               r1                r2                   rout             mix
  ! |   center       |       inside    |    outside          |               |
  ! |<--var_cen      |<--------------- var_in -------------->|     var_out-->|

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: r, z, var_cen, var_in, var_out
    real :: r1, r2, rout, rmix, zmix, zjet, var_zcen, var_zin
    real :: taper
    r1 = sim(nozzle)%rFeatherIn
    r2 = sim(nozzle)%radius
    rout = sim(nozzle)%radius + sim(nozzle)%rFeatherOut
    rmix = sim(nozzle)%radius + sim(nozzle)%rFeatherOut + sim(nozzle)%rFeatherMix
    zjet = sim(nozzle)%length
    zmix = zjet + sim(nozzle)%zFeather

    ! z part
    if (abs(z).ge.0.0 .and. abs(z).lt.zjet) then
      var_zin = 1.0*var_in
      var_zcen = 1.0*var_cen
    else if (abs(z).ge.zjet .and. abs(z).lt.zmix) then
      var_zin = (-2*(zmix-abs(z))**3/(zmix-zjet)**3 &
      + 3*(zmix-abs(z))**2/(zmix-zjet)**2)*(var_in-var_out) + var_out
      var_zcen = (-2*(zmix-abs(z))**3/(zmix-zjet)**3 &
      + 3*(zmix-abs(z))**2/(zmix-zjet)**2)*(var_cen-var_out) + var_out
    else
      var_zin = var_out
      var_zcen = var_out
    end if

    ! radial part
    if (r.ge.0.0 .and. r.lt.r1) then
      taper = (-2*r**3/r1**3 + 3*r**2/r1**2)*(var_zin-var_zcen) + var_zcen
    else if (r.ge.r1 .and. r.lt.rout) then
      taper = 1.0*var_zin
    else if (r.ge.rout .and. r.lt.rmix) then
      taper = (-2*(rmix-r)**3/(rmix-rout)**3 &
      + 3*(rmix-r)**2/(rmix-rout)**2)*(var_zin-var_out) + var_out
    else
      taper = var_out
    endif

  end function taper

  function taperSph(nozzle, r, var_in, var_out)

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: r, var_in, var_out
    real :: rout, rmix
    real :: taperSph
    rout = max(sim(nozzle)%radius+sim(nozzle)%rFeatherOut, sim(nozzle)%length)
    rmix = max( (sim(nozzle)%radius+sim(nozzle)%rFeatherOut+sim(nozzle)%rFeatherMix),&
                (sim(nozzle)%length+sim(nozzle)%zFeather) )

    if (r.ge.0.0 .and. r.lt.rout) then
      taperSph = var_in
    else if (r.ge.rout .and. r.lt.rmix) then
      taperSph = (-2*(rmix-r)**3/(rmix-rout)**3 &
      + 3*(rmix-r)**2/(rmix-rout)**2)*(var_in-var_out) + var_out
    else
      taperSph = var_out
    endif

  end function taperSph

end module Simulation_data


