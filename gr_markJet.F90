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
!!     blk_resolution, which indicates the maximum number of blocks should be
!!     used to resolve one side of the lobes in length. The radial direction
!!     uses half of the blocks, assuming aspect ratio is roughly 2:1.
!!
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
  use Grid_data, ONLY : gr_geometry, gr_smalle
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkIndexLimits
  use Simulation_data, ONLY : sim, sim_onlyHalf
#include "constants.h"
#include "Flash.h"
  implicit none

! Arguments

  integer, intent(IN)   :: nozzle

! Local data

  real, dimension(MDIM) :: blockCenter, blockSize
  integer               :: b, lrefine_0, blk_resolution
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC

  real :: radius, length, sig, distance, theta, vel, fac, pmax, pjet, eintmin
  real, dimension(3) :: plnvec, jetvec, rvec, phivec, velvec
  real, pointer, dimension(:,:,:,:) :: solnDataGC, solnData

  lrefine_0 = min(sim(nozzle)%lrefine_0, lrefine_max)

  ! Maximal number of blocks to resolve one side of the lobes in length
  blk_resolution = 32
  if((gr_geometry == CARTESIAN)) then
     do b = 1, lnblocks
        if(nodetype(b) == LEAF) then
           blockCenter(:) = coord(:,b)
           blockSize(:) = bsize(:,b)

           !call hy_uhd_jetNozzleGeometry(nozzle,blockCenter,radius,length,distance,&
           !                              sig,theta,jetvec,rvec,plnvec,phivec)

           length = sum( sim(nozzle)%coneVec(:)&
                         *( blockCenter(:) - sim(nozzle)%pos(:) ))
           radius = sqrt(sum( (blockCenter(:) - sim(nozzle)%pos(:) &
                               - length*sim(nozzle)%coneVec)**2 ))

           if ( abs(length)/maxval(blockSize) > blk_resolution ) then
              if (lrefine(b) > 1 ) then
                 refine(b)   = .false.
                 derefine(b) = .true.
              else
                 refine(b) = .false.
              endif
           else if ( abs(length)/maxval(blockSize) > blk_resolution/2 ) then
              refine(b) = .false.
           endif
           ! For radial direction, use only half of the blocks, assuming the
           ! aspect ratio is roughly 2:1
           if ( radius/maxval(blockSize) > blk_resolution/2 ) then
              if (lrefine(b) > 1 ) then
                 refine(b)   = .false.
                 derefine(b) = .true.
              else
                 refine(b) = .false.
              endif
           else if ( radius/maxval(blockSize) > blk_resolution/4 ) then
              refine(b) = .false.
           endif
           ! Set the overall maximal refinement level (lrefine_0)
           ! The jet is to be refined by different level (lrefine_max)
           if (lrefine(b) > lrefine_0 ) then
              refine(b)   = .false.
              derefine(b) = .true.
           else if (lrefine(b) == lrefine_0) then
              refine(b) = .false.
           endif
           


           ! Calculate the maximum jet momentum in the block
           call Grid_getBlkIndexLimits(b,blkLimits,blkLimitsGC)
           call Grid_getBlkPtr(b, solnDataGC, CENTER)
           ! Exclude the guard cells!
           solnData => solnDataGC(:, blkLimits(LOW ,IAXIS):blkLimits(HIGH,IAXIS),&
                                     blkLimits(LOW ,JAXIS):blkLimits(HIGH,JAXIS),&
                                     blkLimits(LOW ,KAXIS):blkLimits(HIGH,KAXIS) )
           call Grid_releaseBlkPtr(b, solnDataGC, CENTER)
           pmax = maxval(solnData(JET_SPEC,:,:,:)*solnData(DENS_VAR,:,:,:)*&
                         abs(solnData(VELX_VAR,:,:,:)*sim(nozzle)%coneVec(1)+&
                             solnData(VELY_VAR,:,:,:)*sim(nozzle)%coneVec(2)+&
                             solnData(VELZ_VAR,:,:,:)*sim(nozzle)%coneVec(3)) )
           pjet = sim(nozzle)%velJet*sim(nozzle)%density
           eintmin = minval(solnData(EINT_VAR,:,:,:))

           ! Force maximum refine level for the jet using momentum
           ! Force maximum refinement for abnormally low internal energy region
           ! to increase stability
           if ((pmax >= sim(nozzle)%refine_jetR*pjet) .or. (eintmin <= gr_smalle)) then
              if (lrefine(b) < lrefine_max) then
                 refine(b) = .true.
                 derefine(b) = .false.
              else if (lrefine(b) == lrefine_max) then
                 derefine(b) = .false.
              endif
           endif
           if (pmax >= sim(nozzle)%derefine_jetR*pjet) then
              if (lrefine(b) == lrefine_max) then
                 derefine(b) = .false.
              endif
           endif

           if (sim_onlyHalf) then
              ! Maintain minimal refinement level for lower half of the domain
              if (blockCenter(3) < 0.0) then
                 if (lrefine(b) > 1 ) then
                    refine(b)   = .false.
                    derefine(b) = .true.
                 else
                    refine(b)   = .false.
                 endif
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
