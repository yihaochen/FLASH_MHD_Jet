!!****f* source/physics/sourceTerms/Heat/Heat
!!
!! NAME
!!  
!!  ht_getValueAtPoint
!!
!!
!! SYNOPSIS
!! 
!!  call ht_getValueAtPoint (real(IN), dim(3)             :: pointvec,
!!                             real(OUT), dim(NPROP_VAR)    :: varData)
!!
!!  
!!  
!! DESCRIPTION
!!
!!  Find the value of the variable at a given coordinate.
!!
!!
!! ARGUMENTS
!!
!!  pointvec       : coordinate of the point
!!  varData       : values of the variables of the point
!!
!!***

subroutine ht_getValueAtPoint (blockID, pointvec, del, varData)
!
!==============================================================================
!
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use Simulation_data
  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer,intent(IN) :: blockID
  real, intent(IN), dimension(3) :: pointvec, del
  !real, allocatable, dimension(:,:) :: pointvecs
  real, intent(OUT), dimension(NUNK_VARS) :: varData
  real, pointer, dimension(:,:,:,:) :: solnData
  real, dimension(NUNK_VARS) :: c00,c01,c10,c11,c0,c1

  integer, dimension(1) :: iarr, jarr, karr
  integer :: i0,i1, j0,j1, k0,k1
  real :: xd, yd, zd

  call Grid_getBlkPtr(blockID,solnData,CENTER)
  !allocate(pointvecs(3,size(sim_xCoord)))

  !pointvecs(1,:) = pointvec(1)
  !pointvecs(2,:) = pointvec(2)
  !pointvecs(4,:) = pointvec(3)

  iarr = minloc(abs(sim_xCoord(:)-pointvec(1)))
  jarr = minloc(abs(sim_yCoord(:)-pointvec(2)))
  karr = minloc(abs(sim_zCoord(:)-pointvec(3)))

  !write (*,'(a, 3i11)') 'ijk = ', iarr, jarr, karr
  i0 = iarr(1)
  j0 = jarr(1)
  k0 = karr(1)

  if (sim_xCoord(i0)-pointvec(1) > 0) then
      i0 = i0-1
  endif
  if (sim_yCoord(j0)-pointvec(2) > 0) then
      j0 = j0-1
  endif
  if (sim_zCoord(k0)-pointvec(3) > 0) then
      k0 = k0-1
  endif
  !write (*,'(a, 3es11.3)') 'pnt = ', pointvec
  !write (*,'(a, 3es11.3)') 'xyz0= ', sim_xCoord(i0), sim_yCoord(j0), sim_zCoord(k0)
  i1 = i0+1
  j1 = j0+1
  k1 = k0+1
  !write (*,'(a, 3es11.3)') 'xyz1= ', sim_xCoord(i1), sim_yCoord(j1), sim_zCoord(k1)
  xd = (pointvec(1) - sim_xCoord(i0)) / (del(1))
  yd = (pointvec(2) - sim_yCoord(j0)) / (del(2))
  zd = (pointvec(3) - sim_zCoord(k0)) / (del(3))

  !write (*,'(a, 3es11.3)') 'del = ', del
  !write (*,'(a, 3f11.3)') 'xyzd= ', xd, yd, zd


  ! Interpolation
  c00(:) = solnData(:,i0,j0,k0)*(1.-xd) + solnData(:,i1,j0,k0)*xd
  c10(:) = solnData(:,i0,j1,k0)*(1.-xd) + solnData(:,i1,j1,k0)*xd
  c01(:) = solnData(:,i0,j0,k1)*(1.-xd) + solnData(:,i1,j0,k1)*xd
  c11(:) = solnData(:,i0,j1,k1)*(1.-xd) + solnData(:,i1,j1,k1)*xd

  c0(:) = c00(:)*(1.-yd) + c10*yd
  c1(:) = c01(:)*(1.-yd) + c11*yd

  varData(:) = c0(:)*(1.-zd) + c1(:)*zd
  !write(*,'(3i4)') i0, j0, k0
  !write(*,'(2es11.3)') solnData(DENS_VAR,i0,j0,k0), varData(DENS_VAR)
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)


  return

end subroutine ht_getValueAtPoint


