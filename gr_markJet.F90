!!****if* source/Grid/GridMain/paramesh/gr_markJet
!!
!! NAME
!!  gr_markJet
!!
!!  
!! SYNOPSIS 
!!  gr_markJet(integer(in) :: nozzle) 
!!  
!! PURPOSE 
!!  1. Set the maximal refinement level according to jet heights using parameters
!!  derefine_z1 and derefine_z2 (in unit of the nozzle length). Above each scale 
!!  the max refinement level will decrease by 1. Maximum refinement level below 
!!  derefine_z1 is set by the parameter lrefine_0.
!!  2. Force the jet identified by the momentum to have maximum refinement level.
!!  
!! ARGUMENTS 
!!  nozzle - index of the nozzle (currently 1)
!!  
!! NOTES
!! 
!!  This routine has not yet been tested and should be used only as a guideline for
!!  a user's implementation.
!!  
!!  
!!***

subroutine gr_markJet(nozzle)

!-------------------------------------------------------------------------------
  use tree, ONLY : refine, derefine, lrefine, bsize, coord, lnblocks, nodetype,&
                   lrefine_max
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_geometry
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkIndexLimits
  use Simulation_data
#include "constants.h"
#include "Flash.h"
  implicit none

! Arguments

  integer, intent(IN)   :: nozzle

! Local data

  real, dimension(MDIM) :: blockCenter, blockSize
  integer               :: b, lrefine_0
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC

  real :: radius, length, sig, distance, theta, vel, fac, pmax
  real, dimension(3) :: plnvec, jetvec, rvec, phivec, velvec
  real, pointer, dimension(:,:,:,:) :: solnDataGC, solnData

  lrefine_0 = sim(nozzle)%lrefine_0

  if((gr_geometry == CARTESIAN)) then
     do b = 1, lnblocks
        if(nodetype(b) == LEAF) then
           blockCenter(:) = coord(:,b)
           blockSize(:) = 0.5*bsize(:,b)

           call Grid_getBlkIndexLimits(b,blkLimits,blkLimitsGC)
           call Grid_getBlkPtr(b, solnDataGC, CENTER)
           ! Exclude the guard cells
           solnData => solnDataGC(:, blkLimits(LOW ,IAXIS):blkLimits(HIGH,IAXIS),&
                                     blkLimits(LOW ,JAXIS):blkLimits(HIGH,JAXIS),&
                                     blkLimits(LOW ,KAXIS):blkLimits(HIGH,KAXIS) )
           call Grid_releaseBlkPtr(b, solnDataGC, CENTER)
           pmax = maxval(solnData(JET_SPEC,:,:,:)*solnData(DENS_VAR,:,:,:)*&
                         abs(solnData(VELX_VAR,:,:,:)*jetvec(1)+&
                             solnData(VELY_VAR,:,:,:)*jetvec(2)+&
                             solnData(VELZ_VAR,:,:,:)*jetvec(3)) )

           call hy_uhd_jetNozzleGeometry(nozzle,blockCenter,radius,length,distance,&
                                         sig,theta,jetvec,rvec,plnvec,phivec)

           ! Decrease the maximum refine level when away from the nozzle
           if (abs(length) < sim(nozzle)%derefine_z1*sim(nozzle)%length) then
              if (lrefine(b) > lrefine_0 ) then
                 refine(b)   = .false.
                 derefine(b) = .true.
              else if (lrefine(b) == lrefine_0) then
                 refine(b) = .false.
              endif
           endif
           if (abs(length) >= sim(nozzle)%derefine_z1*sim(nozzle)%length) then
              if (lrefine(b) >= lrefine_0 ) then
                 refine(b)   = .false.
                 derefine(b) = .true.
              else if (lrefine(b) == lrefine_0-1) then
                 refine(b) = .false.
              endif
           endif

           if (abs(length) >= sim(nozzle)%derefine_z2*sim(nozzle)%length) then
              if (lrefine(b) >= lrefine_0-1 ) then
                 refine(b)   = .false.
                 derefine(b) = .true.
              else if (lrefine(b) == lrefine_0-2) then
                 refine(b) = .false.
              endif
           endif
           

           ! Force maximum refine level for the jet using momentum
           if (pmax >= sim(nozzle)%refine_jetR*sim(nozzle)%velocity*sim(nozzle)%density) then
              if (lrefine(b) < lrefine_max) then
                 refine(b) = .true.
                 derefine(b) = .false.
              else if (lrefine(b) == lrefine_max) then
                 derefine(b) = .false.
              endif
           endif
           call Grid_releaseBlkPtr(b, solnData, CENTER)


           
           ! End of leaf-node block loop
        endif
     end do
  else
     call Driver_abortFlash("MarkRefine: geometry other than Cartesian is not yet supported in Jet Simulation")
     !-------------------------------------------------------------------------------
  end if
  return
end subroutine gr_markJet
