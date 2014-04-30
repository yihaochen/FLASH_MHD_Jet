!! NAME
!!
!!  hy_uhd_getA
!!
!! SYNOPSIS
!!
!!  hy_uhd_getA( integer (IN) :: nozzle,
!!               real    (IN) :: simTime,
!!               real    (IN) :: r, 
!!               real    (IN) :: z,
!!               real    (IN) :: phi,
!!               real    (OUT):: Ar,
!!               real    (OUT):: Az,
!!               real    (OUT):: Aphi )
!!
!! DESCRIPTION
!!  
!!  Take the local cylindrical coordinate (r,z,phi) and calculate
!!  the vector potetial at the point. Returns the vector potential
!!  in local cylindrical coordinate.
!!
!!
!! ARGUMENTS
!!
!!  nozzle - index of the nozzle
!!  simTime- simulation time
!!  r      - radial component of the local cylindrical coordinate
!!  z      - z component of the local cylindrical coordinate
!!  phi    - phi angle of the local cylindrical coordinate
!!  Ar     - radial component of the vector potential
!!  Az     - z component of the vector potential
!!  Aphi   - phi component of the vector potential

Subroutine hy_uhd_getA(nozzle,simTime,r,z,phi,Ar,Az,Aphi)
  use Simulation_data

  implicit none

#include "constants.h"
#include "Flash.h"
!#include "UHD.h"

  integer, intent(IN) :: nozzle
  real, intent(IN) :: simTime, r, z, phi
  real :: r1, r2, rjet, vjet
  !real, dimension(3), intent(IN) :: jetvec, rvec, plnvec, phivec 
  real, intent(OUT) :: Ar, Az, Aphi
  !real, dimensiON(3), INTENt(OUT) :: Avec
  INTEGER :: Aopt

  !
  ! geometric factors
  !
  rjet = sim(nozzle)%radius
  r1 = sim(nozzle)%bfeather_inner 
  r2 = rjet - sim(nozzle)%bfeather_outer

  ! velocity of the jet (in z-direction)
  vjet = sim(nozzle)%velocity
  
  !
  ! helical field
  !
  Aopt = 1
  select case(Aopt)
  ! 1) using A_r
    case(1)
      if (r.lt.rjet .and. z.gt.0.0) then
        Ar = (-z + vjet*simTime)*Bphi_model(nozzle, r)
      else if (r.lt.rjet .and. z.lt.0.0) then
        Ar = (-z - vjet*simTime)*Bphi_model(nozzle, r)
      else
        Ar = 0.0
      end if
      Az = 0.0
      Aphi = 0.0

  ! 2) using A_z (divergenless Coulumb gauge) TODO:Time dependence
    case(2)
      if (r.ge.0 .and. r.lt.r1) then
        Az = r**4/(2*r1**3) - r**3/r1**2
      else if (r.ge.r1 .and. r.lt.r2) then
        Az = r
      else if (r.ge.r2 .and. r.lt.rjet) then
        Az = (rjet-r)**4/(2*(rjet-r1)**3) - (rjet-r)**3/(rjet-r1)**2
      else
        Az = 0.0
      end if
      Ar = 0.0
      Aphi = 0.0
    
    case default
        Az = 0.0
        Ar = 0.0
        Aphi = 0.0

  end select

End Subroutine hy_uhd_getA
