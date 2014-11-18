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
!!  Set the maximal refinement level according to jet heights using parameters
!!  derefine_z1 and derefine_z2. Above each scale will decrease max refinement level by 1.
!!  
!! ARGUMENTS 
!!  ic -   Center of the interval/circle/sphere : IAXIS
!!  jc -                                          JAXIS
!!  kc -                                          KAXIS
!!               (Coordinates for nonexistent dimensions are ignored.)
!!  radius -       Radius of the region 
!!  lref  -        If > 0, bring all qualifying blocks to this level of refinement.
!!                 If <= 0, refine qualifying blocks once.
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
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use Simulation_data
#include "constants.h"
#include "Flash.h"
  implicit none

! Arguments

  integer, intent(IN)   :: nozzle

! Local data

  real, dimension(MDIM) :: blockCenter, blockSize
  integer               :: b, lrefine_z0

  real :: radius, length, sig, distance, theta, vel, fac, pmax
  real, dimension(3) :: plnvec, jetvec, rvec, phivec, velvec
  real, pointer, dimension(:,:,:,:) :: solnData

  lrefine_z0 = sim(nozzle)%lrefine_z0

  if((gr_geometry == CARTESIAN)) then
     do b = 1, lnblocks
        if(nodetype(b) == LEAF) then
           blockCenter(:) = coord(:,b)
           blockSize(:) = 0.5*bsize(:,b)

           call hy_uhd_jetNozzleGeometry(nozzle,blockCenter,radius,length,distance,&
                                         sig,theta,jetvec,rvec,plnvec,phivec)

           ! Decrease the maximum refine level when away from the nozzle
           if (abs(length) < sim(nozzle)%derefine_z1*sim(nozzle)%length) then
              if (lrefine(b) > lrefine_z0 ) then
                 refine(b)   = .false.
                 derefine(b) = .true.
              else if (lrefine(b) == lrefine_z0) then
                 refine(b) = .false.
              endif
           endif
           if (abs(length) >= sim(nozzle)%derefine_z1*sim(nozzle)%length) then
              if (lrefine(b) >= lrefine_z0 ) then
                 refine(b)   = .false.
                 derefine(b) = .true.
              else if (lrefine(b) == lrefine_z0-1) then
                 refine(b) = .false.
              endif
           endif

           if (abs(length) >= sim(nozzle)%derefine_z2*sim(nozzle)%length) then
              if (lrefine(b) >= lrefine_z0-1 ) then
                 refine(b)   = .false.
                 derefine(b) = .true.
              else if (lrefine(b) == lrefine_z0-2) then
                 refine(b) = .false.
              endif
           endif
           

           ! Force maximum refine level for the jet using momentum
           call Grid_getBlkPtr(b, solnData, CENTER)
           pmax = maxval(solnData(JET_SPEC,:,:,:)*solnData(DENS_VAR,:,:,:)*&
                         abs(solnData(VELX_VAR,:,:,:)*jetvec(1)+&
                             solnData(VELY_VAR,:,:,:)*jetvec(2)+&
                             solnData(VELZ_VAR,:,:,:)*jetvec(3)) )
           if (pmax >= sim(nozzle)%refine_jetR1*sim(nozzle)%velocity*sim(nozzle)%density) then
              !write (*,'(2i4, 2es11.3)') b, lrefine(b), pmax, sim(nozzle)%velocity*sim(nozzle)%density
              if (lrefine(b) < lrefine_max) then
                 refine(b) = .true.
                 derefine(b) = .false.
              else if (lrefine(b) == lrefine_max) then
                 derefine(b) = .false.
              endif
           endif
           if (pmax >= sim(nozzle)%refine_jetR2*sim(nozzle)%velocity*sim(nozzle)%density) then
              if (lrefine(b) < lrefine_max-1) then
                 refine(b) = .true.
                 derefine(b) = .false.
              else if (lrefine(b) == lrefine_max-1) then
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
