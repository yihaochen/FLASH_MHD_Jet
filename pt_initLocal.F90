!!****if* source/Particles/ParticlesInitialization/WithDensity/CellMassBins/pt_initLocal
!!
!! NAME
!!    pt_initLocal
!!
!! SYNOPSIS
!!
!!    pt_initLocal()
!!
!! DESCRIPTION
!!    Local initialization of  particle locations.  Specific initializations that are
!!      needed only with gas density.  Calculates the total volume and the average density
!!      across all blocks.
!!      Initializes random fields.
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  pt_pRand:                integer  Something to do with initial random distribution
!!
!! NOTES
!!
!!   There is a nice description of the Fortran90 random number routines at
!!         http://www.nsc.liu.se/~boein/f77to90/a5.html#section21c
!!
!!***

!===============================================================================

subroutine pt_initLocal ()

  use Particles_data, ONLY:  pt_geometry,pt_meshMe, pt_meshNumProcs,pt_pRand, &
       pt_totalMass, pt_totalVolume, pt_averageDensity, pt_numParticlesWanted,&
       pt_meshComm

  use Driver_interface, ONLY : Driver_abortFlash

  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr, &
    Grid_getBlkPhysicalSize, Grid_getBlkIndexLimits, &
    Grid_getBlkCenterCoords, Grid_releaseBlkPtr, Grid_getSingleCellVol

  use Simulation_data, ONLY : sim_ptMaxRadius

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"



  ! for random number generator
  integer :: seed_size

  real          :: localMass, localVolume, localDensity
  real          :: dvol
  integer       :: thisBlock, bb, i, j, k, ierr, blockCount, istat
  integer       :: numPEs
  integer       :: blockList(MAXBLOCKS)
  integer, dimension(MDIM) :: point
  integer, dimension(2,MDIM):: blkLimits, blkLimitsGC
  real, dimension(:,:,:,:), pointer :: solnData


  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer :: sizeX,sizeY,sizeZ

  !-------------------------------------------------------------------------------

  ! Runtime Parameters


  ! This incarnation is supposed to support all geometries.


  ! In this routine, we determine the total volume and average density on the
  ! grid and save it.  Note that this will only work correctly if pt_initPositions
  ! has been called after the mesh has been set up OR the DENS_VAR variable contains
  ! accurate zone-average densities.

  localMass = 0.
  localVolume = 0.
  localDensity = 0.

  ! loop over all local leaf blocks
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  do bb = 1, blockCount

     thisBlock = blockList(bb)
     call Grid_getBlkPtr(thisBlock,solnData)

     ! get dimension limits etc.
     call Grid_getBlkIndexLimits(thisBlock,blkLimits,blkLimitsGC)

     sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
     sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
     sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

     allocate(xCoord(sizeX),stat=istat)
     allocate(yCoord(sizeY),stat=istat)
     allocate(zCoord(sizeZ),stat=istat)
     call Grid_getCellCoords(IAXIS,thisBlock,CENTER,.true., xCoord, sizeX)
     call Grid_getCellCoords(JAXIS,thisBlock,CENTER,.true., yCoord, sizeY)
     call Grid_getCellCoords(KAXIS,thisBlock,CENTER,.true., zCoord, sizeZ)


     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        point(3) = k
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           point(2) = j
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              point(1) = i
              call Grid_getSingleCellVol(thisBlock,EXTERIOR,point,dvol)
              if ((xCoord(i)*xCoord(i)+yCoord(j)*yCoord(j)+zCoord(k)*zCoord(k)) &
                  .LE. (sim_ptMaxRadius*sim_ptMaxRadius)) then
                 localMass = localMass + solnData(DENS_VAR,i,j,k)*dvol
                 localVolume = localVolume + dvol
              endif
           enddo
        enddo
     enddo

     !  release the pointer
     call Grid_releaseBlkPtr(thisBlock,solnData)

  enddo  !! of looping over all local leaf blocks
  

  !! get the mass across all processors and store it in data variable pt_totalMass
  call mpi_allreduce(localMass, pt_totalMass, 1, FLASH_REAL, MPI_SUM, &
       pt_meshComm, ierr)
  call mpi_allreduce(localVolume, pt_totalVolume, 1, FLASH_REAL, MPI_SUM, &
       pt_meshComm, ierr)

  ! now calculate the average density
  pt_averageDensity = pt_totalMass / pt_totalVolume

  !-------------------------------------------------------------------------------

  ! randomize the initial particle positions

  !  returned value seed_size gives the number of integers the processor uses for the
  !    starting value
  call random_seed(SIZE=seed_size)
  
  !  generates a large (from pt_pRand) integer, in general different for each processor

  i = int(pt_pRand * pt_meshMe) + pt_meshNumProcs

  !  initializes the random number vector with a fixed seed (from i)
  call random_seed(PUT=(/(i, j = 1, seed_size)/))
  

  !We used to call random_number lots of times but this does not serve
  !any useful purpose for any of our FLASH simulations.


  !-------------------------------------------------------------------------------

  return

end subroutine pt_initLocal


