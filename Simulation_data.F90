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
      real, dimension(3) :: pos, jetvec
      real, dimension(3) :: posOld, jetvecOld
      real, dimension(3) :: angVel, linVel
      real :: radius, length, zTorInj, outflowR
      real :: derefine_z1, derefine_z2, refine_jetR1, refine_jetR2
      integer :: lrefine_0
      real :: power, pressure, density, velocity, gamma, mach
      real :: deltaP, deltaRho
      real :: rFeatherIn, rFeatherOut, zFeather
      real :: beta, helicity
      !TODO: tOn, tOff
      real :: bphi, bz, tOn, tOff, timeMHDon
      real :: bzOld
      real :: t0 = -1.0
      !TODO: initGeometry -> sim_initGeometry
      character(len=MAX_STRING_LENGTH) :: initGeometry

  end TYPE nozzle_struct

  !! *** Runtime Parameters *** !!

  type(nozzle_struct), save, dimension(NOZZLES) :: sim

  real,save :: sim_pAmbient, sim_rhoAmbient, sim_windVel, sim_bzAmbient
  real,save :: sim_densityBeta, sim_densityCoreR
  real,save :: sim_gamma, sim_smallP, sim_smallX
  real,save,allocatable,dimension(:) :: sim_xcoord
  real,save,allocatable,dimension(:) :: sim_ycoord
  real,save,allocatable,dimension(:) :: sim_zcoord
  real,save,allocatable,dimension(:) :: sim_xcoordf
  real,save,allocatable,dimension(:) :: sim_ycoordf
  real,save,allocatable,dimension(:) :: sim_zcoordf
  character(len=MAX_STRING_LENGTH) :: sim_densityProfile


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
    r1 = sim(nozzle)%rFeatherIn
    r2 = sim(nozzle)%radius
    rout = sim(nozzle)%radius + sim(nozzle)%rFeatherOut

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

  function taperL(nozzle, z, var_in, var_out)

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: z, var_in, var_out
    real :: zout, zjet
    real :: taperL
    zjet = sim(nozzle)%length
    zout = zjet + sim(nozzle)%zFeather

    ! z part
    if (abs(z).ge.0.0 .and. abs(z).lt.zjet) then
      taperL = 1.0*var_in
    else if (abs(z).ge.zjet .and. abs(z).lt.zout) then
      taperL = (-2*(zout-abs(z))**3/(zout-zjet)**3 &
      + 3*(zout-abs(z))**2/(zout-zjet)**2)*(var_in-var_out) + var_out
    else
      taperL = var_out
    end if

  end function taperL
  
  function taper(nozzle, r, z, var_cen, var_in, var_out)
  ! Taper function for both R and z direction
  ! For radial direction, variable could have 3 values: center, inside, and outside
  !
  ! 0          rFeatherIn       radius         radius+rFeatherOut
  ! 0               r1                r2                   rout
  ! |   center       |       inside    |    outside          |
  ! |   var_cen      |       var_in    |    var_out          |

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: r, z, var_cen, var_in, var_out
    real :: r1, r2, rout, zout, zjet, var_zcen, var_zin
    real :: taper
    r1 = sim(nozzle)%rFeatherIn
    r2 = sim(nozzle)%radius
    rout = sim(nozzle)%radius + sim(nozzle)%rFeatherOut
    zjet = sim(nozzle)%length
    zout = zjet + sim(nozzle)%zFeather

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

  function taperSph(nozzle, r, var_in, var_out)

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: r, var_in, var_out
    real :: r2, rout
    real :: taperSph
    r2 = max(sim(nozzle)%radius, sim(nozzle)%length)
    rout = max( (sim(nozzle)%radius+sim(nozzle)%rFeatherOut),&
                (sim(nozzle)%length+sim(nozzle)%zFeather) )

    if (r.ge.0.0 .and. r.lt.r2) then
      taperSph = var_in
    else if (r.ge.r2 .and. r.lt.rout) then
      taperSph = (-2*(rout-r)**3/(rout-r2)**3 &
      + 3*(rout-r)**2/(rout-r2)**2)*(var_in-var_out) + var_out
    else
      taperSph = var_out
    endif

  end function taperSph

end module Simulation_data


