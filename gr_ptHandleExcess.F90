!!****if* source/Grid/GridParticles/gr_ptHandleExcess
!!
!! NAME
!!
!!  gr_ptHandleExcess
!!
!! SYNOPSIS
!!
!!  gr_ptHandleExcess(real,dimension(propCount),intent(INOUT) :: particles,
!!                            integer,intent(IN)    :: propCount,
!!                            integer,intent(INOUT) :: localNum,
!!                            integer,intent(IN) :: maxPerProc)
!!
!! DESCRIPTION
!!   This routine is called if the runtime parameter gr_ptRemove
!!   is activated. This particular implementation is meant to be
!!   mostly an example for users, since they are likely to want their
!!   own algorithm for removal of particles.
!!
!! ARGUMENTS
!!
!!   particles : Data structure holding particles
!!
!!   propCount : the count of fields in the particles data structure
!!
!!   localNum  : number of particles on a processor
!!
!!   maxPerProc : maximum number of particles allowed on a processor
!!
!!
!!***

!!REORDER(4): solnData

#ifdef DEBUG_ALL
#define DEBUG_GRIDPARTICLES
#endif

subroutine gr_ptHandleExcess(particles,propCount,localNum,maxPerProc)
  use gr_ptData, ONLY : gr_ptRemove,gr_ptRemoveAlgo,gr_ptNumToReduce,&
       gr_ptBlk, gr_ptTag
  use Grid_interface, ONLY : Grid_getListOfBlocks,Grid_getBlkPtr,Grid_releaseBlkPtr,&
       Grid_getBlkIndexLimits
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

#include "Flash.h"
#include "constants.h"

  integer,intent(IN) :: propCount
  integer,intent(IN) :: maxPerProc
  integer,intent(INOUT) :: localNum
  real,dimension(propCount,maxPerProc),intent(INOUT) :: particles
  real,dimension(MAXBLOCKS) :: rho
  integer,dimension(MAXBLOCKS) :: blkList
  integer :: i,j,k,prop,blockID,blkCount,maxRhoBlock
  integer,dimension(MDIM) :: ind
  real :: maxRho, minRho, tol
  integer,dimension(LOW:HIGH,MDIM)::blkLimits,blkLimitsGC

  integer :: localNumIn, numReduced
  integer,save :: digit = 10
  real,dimension(:,:,:,:),pointer :: solnData

  if(gr_ptRemove) then
#ifdef DEBUG_GRIDPARTICLES
     print*,'called gr_ptHandleExcess(particles,',propCount,localNum,maxPerProc,')'
#endif
     select case(gr_ptRemoveAlgo)
     case(1) !! remove particles based upon some density values
#ifdef DENS_VAR        
        call Grid_getListOfBlocks(LEAF,blkList,blkCount)
        tol=0.9
        maxRho=0.0
        ind=1
        do i = 1,blkCount
           blockID=blkList(i)
           call Grid_getBlkPtr(blockID,solnData,CENTER)
           call Grid_getBlkIndexLimits(blockID,blkLimitsGC,blkLimits)
           ind(1:NDIM)=blkLimits(LOW,1:NDIM)-blkLimitsGC(LOW,1:NDIM)+&
                (blkLimits(HIGH,1:NDIM)- blkLimits(LOW,1:NDIM)+1)/2
           rho(i)=solnData(DENS_VAR,ind(IAXIS),ind(JAXIS),ind(KAXIS))
           if(maxRho<rho(i)) then
              maxRho=rho(i)
              maxRhoBlock=blockID
           end if
           call Grid_releaseBlkPtr(blockID,solnData)
        end do
        
        minRho=minval(rho(1:blkCount))
        
        if((minRho/maxRho)>tol) then
           
           blockID=1
           do while (blockID<blkCount)
              k=localNum
              i=1
              do j = 1,localNum
                 prop = int(particles(gr_ptBlk,i))
                 if(prop==blkList(blockID)) then
                    particles(1:propCount,i)=particles(1:propCount,k)
                    k=k-1
                    blockID=blockID+1
                    
                 else
                    i=i+1
                 end if
              end do
              localNum=i-1
              
              blockID=blockID+1
           end do
        else
           
           numReduced=0
           i=1
           k=localNum
           
           do while((numReduced < gr_ptNumToReduce).and.(i<localNum))
              prop=int(particles(gr_ptBlk,i))
              if(prop==maxRhoblock) then
                 particles(1:propCount,i)=particles(1:propCount,k)
                 k=k-1
                 numReduced=numReduced+1
              else
                 i=i+1
              end if
           end do
           localNum=i-1
        end if
#endif        
     case(2) !! remove particles with certain digits in their tag number
        localNumIn = localNum;  i = localNum+1 
        numReduced = 0
        do while ((localNum > maxPerProc .OR. numReduced < gr_ptNumToReduce) &
                  .AND. (localNum > 0))
           digit=mod(digit+9,10)   !9,8,7,6,5,4,3,2,1,0,9,8,7,6,..
#ifdef DEBUG_GRIDPARTICLES
           print*,'tags',particles(gr_ptTag,1:localNum),gr_ptTag,digit
#endif
           k=localNum
           i=1
           do j=1,localNum
              prop = int(particles(gr_ptTag,i))
              if(mod(prop,10)==digit) then
                 particles(1:propCount,i)=particles(1:propCount,k)
                 k=k-1
                 numReduced=numReduced+1
              else
                 i=i+1
              end if
           end do
           print*,'started with ',localNum,' now have ',i-1,'...'
           localNum=i-1
        end do
        print*,'started with ',localNumIn,' now have ',i-1,'.'
     end select
  else
     call Driver_abortFlash("The number of particles exceeded maximum allowed")
  end if
  return
end subroutine gr_ptHandleExcess
