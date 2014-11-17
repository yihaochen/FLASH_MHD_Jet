!!****if* source/Grid/GridMain/paramesh/gr_markJetHeight
!!
!! NAME
!!  gr_markJetHeight
!!
!!  
!! SYNOPSIS 
!!  gr_markJetHeight(integer(in) :: nozzle) 
!!  
!! PURPOSE 
!!  Set the maximal refinement level according to jet heights using parameters
!!  deref_z1 and deref_z2. Above each scale will decrease max refinement level by 1.
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

subroutine gr_markJetHeight(nozzle)

!-------------------------------------------------------------------------------
  use tree, ONLY : refine, derefine, lrefine, bsize, coord, lnblocks, nodetype,\
                   lrefine_max
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_geometry
  use Simulation_data
#include "constants.h"
#include "Flash.h"
  implicit none

! Arguments

  integer, intent(IN)   :: nozzle

! Local data

  real, dimension(MDIM) :: blockCenter, blockSize
  integer               :: b

  real :: radius, length, sig, distance, theta, vel, fac
  real, dimension(3) :: plnvec, jetvec, rvec, phivec, velvec


  if((gr_geometry == CARTESIAN)) then
     do b = 1, lnblocks
        if(nodetype(b) == LEAF) then
           blockCenter(:) = coord(:,b)
           blockSize(:) = 0.5*bsize(:,b)

           call hy_uhd_jetNozzleGeometry(nozzle,blockCenter,radius,length,distance,&
                                         sig,theta,jetvec,rvec,plnvec,phivec)

           ! Force maximum refine level within the nozzle 
           if ((abs(length) <= sim(nozzle)%length) .and. (radius <= sim(nozzle)%radius)) then
              if (lrefine(b) < lrefine_max) then
                 refine(b) = .true.
                 derefine(b) = .false.
              else if (lrefine(b) == lrefine_max) then
                 derefine(b) = .false.
              endif
           endif
           
           ! Decrease the maximum refine level when away from the nozzle
           if (abs(length) >= sim(nozzle)%deref_z1*sim(nozzle)%length) then
              if (lrefine(b) >= lrefine_max ) then
                 refine(b)   = .false.
                 derefine(b) = .true.
              else if (lrefine(b) == lrefine_max-1) then
                 refine(b) = .false.
              endif
           endif

           if (abs(length) >= sim(nozzle)%deref_z2*sim(nozzle)%length) then
              if (lrefine(b) >= lrefine_max-1 ) then
                 refine(b)   = .false.
                 derefine(b) = .true.
              else if (lrefine(b) == lrefine_max-2) then
                 refine(b) = .false.
              endif
           endif
           
           ! End of leaf-node block loop
        endif
     end do
  else
     call Driver_abortFlash("MarkRefine: geometry other than Cartesian is not yet supported in Jet Simulation")
     !-------------------------------------------------------------------------------
  end if
  return
end subroutine gr_markJetHeight
