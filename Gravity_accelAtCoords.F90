!!****if* source/physics/Gravity/GravityMain/PointMass/Gravity_accelAtCoords
!!
!! NAME
!!
!!  Gravity_accelAtCoords 
!!
!! SYNOPSIS
!!
!!  Gravity_accelAtCoords(integer(IN) :: numPoints,
!!                      real(IN)      :: iCoords(:),
!!                      real(IN)      :: jCoords(:),
!!                      real(IN)      :: kCoords(:),
!!                      integer(IN)   :: accelDir,
!!                      real(OUT)     :: accel(numPoints),
!!                      integer(IN)   :: blockID,
!!                      integer(IN),optional :: potentialIndex)
!!
!! DESCRIPTION
!!
!!  This routine computes the gravitational acceleration in a
!!  specified direction for a vector of points given by their
!!  coordinates.
!!
!! ARGUMENTS
!!
!!  iCoords,jCoords,kCoords: coordinates of the points where the
!!                           gravitational accelation is requested.
!!                           Each of these arrays should either be
!!                           of lenght numPoints (or more), in which
!!                           case its nth value is used for the nth
!!                           point; or else of dimension 1, in which
!!                           case the value is used for all points.
!!  accelDir :    The acceleration direction:  allowed values are 
!!              IAXIS, JAXIS and IAXIS. These values are defined
!!              in constants.h.
!!  numPoints :  Number of cells to update in accel()
!!  accel     :   Array to receive results
!!  blockID  :  The local identifier of the block to work on,
!!                not applicable in pointmass gravity.
!!  potentialIndex :  optional, not applicable in pointmass gravity
!! 
!!***

subroutine Gravity_accelAtCoords (numPoints, iCoords,jCoords,kCoords, accelDir,&
     accel, blockID, &
     potentialIndex)

!=======================================================================

  use Gravity_data, ONLY: grv_ptxpos, grv_ptypos, grv_ptzpos, grv_factor, &
       useGravity
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use Simulation_data

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: accelDir, numPoints
  real, dimension(:),INTENT(in) :: iCoords,jCoords,kCoords
  real, dimension(numPoints),INTENT(OUT) :: accel
  integer, intent(IN),optional :: blockID
  integer, intent(IN),optional :: potentialIndex

!==========================================================================

#ifdef FIXEDBLOCKSIZE
  real,dimension(numPoints) ::xCenter,yCenter,zCenter
#else
  real,allocatable,dimension(:) ::xCenter,yCenter,zCenter
#endif
  real :: r, gasConst

  integer :: ii, nozzle=1

!==============================================================================
  if (.NOT.useGravity .or. sim_densityProfile=="uniform") then
     accel(1:numPoints) = 0.0
     return
  end if

  call PhysicalConstants_get("ideal gas constant", gasConst)
#ifndef FIXEDBLOCKSIZE
  allocate(xCenter(numPoints))
  allocate(yCenter(numPoints))
  allocate(zCenter(numPoints))
#endif
  zCenter = 0.
  yCenter = 0.
  if (NDIM == 3) then 
     if (size(kCoords) .GE. numPoints) then
        zCenter(1:numPoints) = kCoords(1:numPoints) - sim(nozzle)%pos(3)
     else
        zCenter(1:numPoints) = kCoords(1) - sim(nozzle)%pos(3)
     end if

  endif
  if (NDIM >= 2) then
     if (size(jCoords) .GE. numPoints) then
        yCenter(1:numPoints) = jCoords(1:numPoints) - sim(nozzle)%pos(2)
     else
        yCenter(1:numPoints) = jCoords(1) - sim(nozzle)%pos(2)
     end if

  endif
  if (size(iCoords) .GE. numPoints) then
     xCenter = iCoords(1:numPoints) - sim(nozzle)%pos(1)
  else
     xCenter = iCoords(1) - sim(nozzle)%pos(1)
  end if

  if (accelDir .eq. IAXIS) then                       ! x-component
     do ii = 1, numPoints
        r=sqrt(xCenter(ii)*xCenter(ii) + yCenter(ii)*yCenter(ii) + zCenter(ii)*zCenter(ii))
        accel(ii) = gasConst/sim_mu*xCenter(ii)&
                    ! Density derivative
                    *(-3.*sim_densityBeta/(1.+r*r/sim_rCore**2)/sim_rCore**2&
                    ! Temperature profile
                    *sim_Tout*(1.0+(r/sim_rCoreT)**3)&
                    /(sim_Tout/sim_Tcore+(r/sim_rCoreT)**3)&
                    ! Temperature derivative
                    +3.0*r*sim_Tout*(sim_Tout/sim_Tcore-1.0)*(sim_rCoreT)**3&
                    /(sim_Tout/sim_Tcore*(sim_rCoreT)**3+r**3)**2)

     end do

  else if (accelDir .eq. JAXIS) then          ! y-component

     do ii = 1, numPoints
        r=sqrt(xCenter(ii)*xCenter(ii) + yCenter(ii)*yCenter(ii) + zCenter(ii)*zCenter(ii))
        accel(ii) = gasConst/sim_mu*yCenter(ii)&
                    ! Density derivative
                    *(-3.*sim_densityBeta/(1.+r*r/sim_rCore**2)/sim_rCore**2&
                    ! Temperature profile
                    *sim_Tout*(1.0+(r/sim_rCoreT)**3)&
                    /(sim_Tout/sim_Tcore+(r/sim_rCoreT)**3)&
                    ! Temperature derivative
                    +3.0*r*sim_Tout*(sim_Tout/sim_Tcore-1.0)*(sim_rCoreT)**3&
                    /(sim_Tout/sim_Tcore*(sim_rCoreT)**3+r**3)**2)
     end do

  else if (accelDir .eq. KAXIS) then          ! z-component

     do ii = 1, numPoints
        r=sqrt(xCenter(ii)*xCenter(ii) + yCenter(ii)*yCenter(ii) + zCenter(ii)*zCenter(ii))
        accel(ii) = gasConst/sim_mu*zCenter(ii)&
                    ! Density derivative
                    *(-3.*sim_densityBeta/(1.+r*r/sim_rCore**2)/sim_rCore**2&
                    ! Temperature profile
                    *sim_Tout*(1.0+(r/sim_rCoreT)**3)&
                    /(sim_Tout/sim_Tcore+(r/sim_rCoreT)**3)&
                    ! Temperature derivative
                    +3.0*r*sim_Tout*(sim_Tout/sim_Tcore-1.0)*(sim_rCoreT)**3&
                    /(sim_Tout/sim_Tcore*(sim_rCoreT)**3+r**3)**2)
     end do

  end if

!==============================================================================
#ifndef FIXEDBLOCKSIZE
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
#endif

  return

end subroutine Gravity_accelAtCoords
