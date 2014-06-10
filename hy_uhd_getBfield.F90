!! NAME
!!
!!  hy_uhd_getBfield
!!
!! SYNOPSIS
!!
!!  hy_uhd_getBfield( integer (IN) :: nozzle,
!!                    real    (IN) :: simTime,
!!                    real    (IN) :: r, 
!!                    real    (IN) :: z,
!!                    real    (IN) :: phi,
!!                    real    (OUT):: Br,
!!                    real    (OUT):: Bz,
!!                    real    (OUT):: Bphi )
!!
!! DESCRIPTION
!!  
!!  Take the local cylindrical coordinate (r,z,phi) and calculate
!!  the magnetic field from the vector potential (hy_uhd_getA.F90)
!!  at a given point. Returns the magnetic field  in local 
!!  cylindrical coordinate.
!!
!!
!! ARGUMENTS
!!
!!  nozzle - index of the nozzle
!!  simTime- simulation time
!!  r      - radial component of the local cylindrical coordinate
!!  z      - z component of the local cylindrical coordinate
!!  phi    - phi angle of the local cylindrical coordinate
!!  Br     - radial component of the magnetic field
!!  Bz     - z component of the magnetic field
!!  Bphi   - phi component of the magnetic field

Subroutine hy_uhd_getBfield(nozzle,simTime,r,z,phi,Br,Bz,Bphi)
  use Simulation_data

  implicit none

#include "constants.h"
#include "Flash.h"
!#include "UHD.h"

  integer, intent(IN) :: nozzle
  real, intent(IN) :: simTime, r, z, phi
  real, intent(OUT) :: Br, Bz, Bphi
  real :: Ar1, Ar2, Az, Aphi, del


  Bphi = taper(nozzle, r, sim(nozzle)%bphi, 0.0)
  Br = 0.0
  Bz = 0.0

End Subroutine hy_uhd_getBfield
