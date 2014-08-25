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
  real :: r1, r2, rout, rpol, bf, vjet, const
  !real, dimension(3), intent(IN) :: jetvec, rvec, plnvec, phivec 
  real, intent(OUT) :: Ar, Az, Aphi
  !real, dimensiON(3), INTENt(OUT) :: Avec
  INTEGER :: Aopt

  ! for convenience
  bf = sim(nozzle)%bfeather_outer
  !
  ! geometric factors
  !
  r1 = sim(nozzle)%bfeather_inner 
  r2 = sim(nozzle)%radius
  rout = sim(nozzle)%radius + bf
  rpol = rout + bf


  ! velocity of the jet (in z-direction)
  vjet = sim(nozzle)%velocity
  
  !
  ! toroidal field
  !
  Aopt = 2
  select case(Aopt)
  ! 1) using A_r
    case(1)
      if (r.lt.rout) then
        !Ar = (-z + vjet*simTime)*taperR(nozzle, r, sim(nozzle)%bphi, 0.0)
        Ar = coshat(z, sim(nozzle)%bPosZ, sim(nozzle)%zfeather)*&
             taperR(nozzle, r, 1.0, 0.0)*sim(nozzle)%bphi
      else
        Ar = 0.0
      end if
      Az = 0.0
      Aphi = 0.5*r*sim(nozzle)%bz

  ! 2) using A_z (divergenless Coulumb gauge?) 
    case(2)
      const = r2*(0.5*r2-4.0/PI**2/r2*bf**2)
      if (r.ge.0 .and. r.le.r1) then
        Az = r**4/(2*r1**3) - r**3/r1**2 + 0.5*(-r1+rout+r2)
        Aphi = 0.5*r
      else if (r.gt.r1 .and. r.le.r2) then
        Az = -r + 0.5*(rout+r2)
        Aphi = 0.5*r
      else if (r.gt.r2 .and. r.le.rout) then
        Az = -(rout-r)**4/(2*(rout-r2)**3) + (rout-r)**3/(rout-r2)**2
        !Aphi = 0.5*sim(nozzle)%bz*r2
        Aphi = (2.0/PI**2*bf*(2*bf/r*sin(PI/2.0/bf*(rout-r))+PI*cos(PI/2.0/bf*(rout-r))) + const/r)*&
                0.5*(1.0+cos(PI*(max(rout,min(rpol,r))-rout)/(rpol-rout)))
      else if (r.gt.rout .and. r.le.rpol) then
        Az = 0.0
        !Aphi = 0.5*sim(nozzle)%bz*r2
        Aphi = (2.0/PI**2*bf*(2*bf/r*sin(PI/2.0/bf*(rout-r))+PI*cos(PI/2.0/bf*(rout-r))) + const/r)*&
                0.5*(1.0+cos(PI*(max(rout,min(rpol,r))-rout)/(rpol-rout)))
      else
        Az = 0.0
        Aphi = 0.0
      end if

      Az = Az*sim(nozzle)%bphi*sim(nozzle)%velocity/sim(nozzle)%zfeather*&
           0.5*(1.0+cos(PI*max(-1.0,(min(1.0,(abs(z)-sim(nozzle)%bPosZ)/sim(nozzle)%zfeather)))))

      Ar = 0.0

      Aphi = Aphi*sim(nozzle)%bz*&
           0.5*(1.0+cos(PI*max(0.0,(min(1.0,(abs(z)-sim(nozzle)%bPosZ)/sim(nozzle)%zfeather)))))
      ! exponential decay of Aphi results in the discontinuity of B_r
      !Aphi = Aphi*exp(-(max(0.0, abs(z)-sim(nozzle)%bPosZ))/sim(nozzle)%zfeather)
    
    case default
      Az = 0.0
      Ar = 0.0
      Aphi = 0.0

  end select

End Subroutine hy_uhd_getA
