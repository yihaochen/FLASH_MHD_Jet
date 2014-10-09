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
      real :: radius, length, bPosZ
      real :: power, pressure, density, velocity, gamma, mach
      real :: deltaP, deltaRho
      real :: rfeather_inner, rfeather_outer, zfeather
      real :: bphi, bz, timeMHDon
  end TYPE nozzle_struct

  !! *** Runtime Parameters *** !!

  type(nozzle_struct), save, dimension(NOZZLES) :: sim

  real,save :: t0 = -1.0
  real,save :: sim_pAmbient, sim_rhoAmbient, sim_windVel, sim_bzAmbient
  real,save :: sim_gamma, sim_smallP, sim_smallX
  real,save,allocatable,dimension(:) :: sim_xcoord
  real,save,allocatable,dimension(:) :: sim_ycoord
  real,save,allocatable,dimension(:) :: sim_zcoord
  real,save,allocatable,dimension(:) :: sim_xcoordf
  real,save,allocatable,dimension(:) :: sim_ycoordf
  real,save,allocatable,dimension(:) :: sim_zcoordf


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
    r1 = sim(nozzle)%rfeather_inner
    r2 = sim(nozzle)%radius
    rout = sim(nozzle)%radius + sim(nozzle)%rfeather_outer

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
    zout = zjet + sim(nozzle)%zfeather

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
  ! 0          rfeather_inner       radius         radius+rfeather_outer
  ! 0               r1                r2                   rout
  ! |   center       |       inside    |    outside          |
  ! |   var_cen      |       var_in    |    var_out          |

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: r, z, var_cen, var_in, var_out
    real :: r1, r2, rout, zout, zjet, var_zcen, var_zin
    real :: taper
    r1 = sim(nozzle)%rfeather_inner
    r2 = sim(nozzle)%radius
    rout = sim(nozzle)%radius + sim(nozzle)%rfeather_outer
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


  function ETaper(nozzle, r, z, var_cen, var_in, var_out)
  ! Taper function for electric field in both R and z direction
  ! The electric nozzle is bigger than the hydro nozzle to avoid boundary
  ! magnification of the magnetic field due to the large velocity.
  ! See taper() for more details.

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: r, z, var_cen, var_in, var_out
    real :: r1, r2, rout, zout, zjet, var_zcen, var_zin
    real :: ETaper
    r1 = sim(nozzle)%rfeather_inner
    r2 = sim(nozzle)%radius + sim(nozzle)%rfeather_outer
    rout = r2 + sim(nozzle)%rfeather_outer
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
      ETaper = (-2*r**3/r1**3 + 3*r**2/r1**2)*(var_zin-var_zcen) + var_zcen
    else if (r.ge.r1 .and. r.lt.r2) then
      ETaper = 1.0*var_zin
    else if (r.ge.r2 .and. r.lt.rout) then
      ETaper = (-2*(rout-r)**3/(rout-r2)**3 &
      + 3*(rout-r)**2/(rout-r2)**2)*(var_zin-var_out) + var_out
    else
      ETaper = var_out
    endif

  end function ETaper


  function coshat(x, halfwidth, feather, max_in)
  !                _____________________________ 
  !              /                               \
  !          /                                       \
  !      /                                               \
  ! ---------------------------------------------------------
  !    |          |              |              |          |
  !          |     halfwidth     |
  !    |  feather  |

    real, INTENT(in) :: x, halfwidth, feather
    real, INTENT(in), optional :: max_in
    real :: coshat, hatvalue

    if (present(max_in)) then
        hatvalue = max_in
    else
        hatvalue = 2.0*feather/PI
    endif
    
    if (abs(x).le.halfwidth-0.5*feather) then
        coshat = hatvalue
    else if (abs(x).lt.(halfwidth+0.5*feather)) then
        coshat = hatvalue*0.5*(1.0+cos( PI*((abs(x) - halfwidth)/feather+0.5) ))
    else
        coshat = 0.0
    endif

  end function coshat

  subroutine calc_jet(nozzle, time, dt)
    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: time, dt
    ! Calculate the jet pressure
    real :: p, g, v, r, L
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

  end subroutine calc_jet

end module Simulation_data


