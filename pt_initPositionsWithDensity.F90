!!****if* source/Particles/ParticlesInitialization/WithDensity/CellMassBins/pt_initPositionsWithDensity
!!
!! NAME
!!
!!  pt_initPositionsWithDensity
!!
!! SYNOPSIS
!!
!!  call pt_initPositionsWithDensity(integer(IN)  :: blockID,
!!                                   logical(OUT) :: success)
!!
!! DESCRIPTION
!!
!!    Initialize particle locations.  This version sets up passive tracer
!!      particles that are distributed randomly according to the gas density     
!!
!! ARGUMENTS
!!
!!   blockID:  ID of block in current processor
!!
!!   success:       returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!! PARAMETERS
!!
!!  pt_numParticlesWanted:   integer  Approximate number of tracer particles to use
!!                                throughout domain ??
!!
!!***
!========================================================================

subroutine pt_initPositionsWithDensity (blockID,success)

  use Particles_data, ONLY:  pt_meshMe, pt_numLocal, &
       pt_maxPerProc, particles, pt_velNumAttrib, pt_velAttrib,pt_typeInfo
  use Particles_data, ONLY : pt_numParticlesWanted, pt_totalMass

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_getSingleCellVol, Grid_getDeltas, &
    Grid_getBlkBoundBox, Grid_getCellCoords
  use Driver_interface, ONLY:  Driver_abortFlash
  use Logfile_interface, ONLY:  Logfile_stamp
  use Particles_interface, ONLY : Particles_mapFromMesh
  use pt_interface, ONLY : pt_initLocal
  use Simulation_data, ONLY : sim_ptMaxRadius
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"
  integer, INTENT(in) :: blockID
  logical, INTENT(out) :: success

  integer          :: newParticlesThisBlock
  real          :: blockMass
  real          :: dvol
  real          :: xpos, ypos, zpos, xvel, yvel, zvel, radius
  integer       :: p, b, bb, i, j, k, numParticlesThisBlock, tag
  integer       :: pAdd

  real, dimension(:,:,:,:), pointer :: solnData

  real    :: urp  ! random number generator
  integer :: part_props=NPART_PROPS
  integer :: myType=1
  integer :: partID, ind
  integer :: blockList(MAXBLOCKS)
  integer :: iSize,jSize,kSize,iSizeGC,jSizeGC,kSizeGC,numCells,icell
  real :: accu
  integer, dimension(MDIM) :: point  ! indices for a given location
  real, dimension(MDIM) :: delta, pos,partAttribVec
  integer, dimension(2,MDIM):: blkLimits, blkLimitsGC
  real,dimension(:),allocatable :: xLeftCoord,yLeftCoord,zLeftCoord
  real,dimension(:,:,:),allocatable :: cellVol
  integer,dimension(:,:),allocatable :: icell2ijk
  real,dimension(:),allocatable :: mass, cumuMass
  real,dimension(LOW:HIGH,MDIM) :: bndBox

  !----------------------------------------------------------------

  ! Access the mesh data for this block.
  call Grid_getBlkPtr(blockID,solnData)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getDeltas(blockID,delta)
  call Grid_getBlkBoundBox(blockID,bndBox)
  iSize = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
  jSize = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
  kSize = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1
  numCells = iSize*jSize*kSize
  iSizeGC = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  jSizeGC = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  kSizeGC = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(xLeftCoord(iSizeGC))
  allocate(yLeftCoord(jSizeGC))
  allocate(zLeftCoord(kSizeGC))
  allocate(cellVol(iSizeGC,jSizeGC,kSizeGC))

  !----------------------------------------------------------------------
  call Grid_getCellCoords(IAXIS, blockID, LEFT_EDGE, .true., xLeftCoord, iSizeGC)
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS, blockID, LEFT_EDGE, .true., yLeftCoord, jSizeGC)
  end if
  if (NDIM == 3) then
     call Grid_getCellCoords(KAXIS, blockID, LEFT_EDGE, .true., zLeftCoord, kSizeGC)
  end if

!  call Grid_getBlkData(blockID, CELL_VOLUME, -1, EXTERIOR, (/1, 1, 1/), cellVol, (/iSize, jSize, kSize/))

  ! Compute the amount of mass in this block and find the maximum density.

  blockMass = 0.

! This duplicates the loop below, and could perhaps be consolidated with it.
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     point(3) = k
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        point(2) = j
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           point(1) = i
           call Grid_getSingleCellVol(blockID,EXTERIOR,point,dvol)
           cellVol(i,j,k) = dvol
           blockMass = blockMass + solnData(DENS_VAR,i,j,k)*dvol
        enddo
     enddo
  enddo

  numParticlesThisBlock = int(anint(pt_numParticlesWanted * (blockMass / pt_totalMass)))

  !  Check that this requested number isn't going to blow data allocation
  success=.true.
  if (pt_numLocal + numParticlesThisBlock > pt_maxPerProc) then
     call Logfile_stamp(numParticlesThisBlock,"This block mass would generate too many additional particles:")
     call Grid_releaseBlkPtr(blockID,solnData) !cleanup before premature return...
     deallocate(cellVol)
     deallocate(zLeftCoord)
     deallocate(yLeftCoord)
     deallocate(xLeftCoord)
     success=.false.
     pt_numLocal=0
     return
  end if
  
  allocate(icell2ijk(MDIM,numCells))
  allocate(mass(numCells))
  allocate(cumuMass(numCells))


  ! Initialize the particles for this block.



  mass = 0.0
  accu = 0.0
  icell = 0
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           icell = icell+1
           icell2ijk(IAXIS,icell) = i
           icell2ijk(JAXIS,icell) = j
           icell2ijk(KAXIS,icell) = k
           mass(icell) = solnData(DENS_VAR,i,j,k) * cellVol(i,j,k)
           accu = accu + mass(icell)
           cumuMass(icell) = accu
        end do
     end do
  end do
  if (icell .NE. numCells) then
     call Driver_abortFlash('pt_initPositionsWithDensity: Garbage icell after accumulation loop!')
  end if
  deallocate(cellVol)

  pAdd = 0
  do p = 1, numParticlesThisBlock


     !  NOTE that some whole complicated (and LBR doesn't understand it) procedure
     !  has been called in Particles_initWithDensity to initialize the random seed.

     call random_number(urp)
     urp = urp * blockMass
     call ut_hunt(cumuMass,numCells,urp,icell)
     icell = icell + 1
     if (icell > numCells) then
        call Driver_abortFlash('pt_initPositionsWithDensity: Garbage after ut_hunt call!')
     end if


     i = icell2ijk(IAXIS,icell)
     j = icell2ijk(JAXIS,icell)
     k = icell2ijk(KAXIS,icell)

     ! Calculate X position within cell, needed in all geometries
     call random_number (harvest=urp)

     xpos = urp * delta(IAXIS)             
     xpos = xpos + xLeftCoord(i)

     if (NDIM >= 2) then
        call random_number (harvest=urp)
        ypos = urp * delta(JAXIS)
        ypos = ypos + yLeftCoord(j)
     else
        ypos = 0.
     endif

     if (NDIM == 3) then
        call random_number (harvest=urp)
        zpos = urp * delta(KAXIS)
        zpos = zpos + zLeftCoord(k)
     else
        zpos = 0.
     endif

     if ((xpos*xpos+ypos*ypos+zpos*zpos) .LE. (sim_ptMaxRadius*sim_ptMaxRadius)) then

        pAdd = pAdd+1

        pos(IAXIS)=xpos
        pos(JAXIS)=ypos
        pos(KAXIS)=zpos

           
        ! now set up 
        partID=pt_numLocal+pAdd

        particles(BLK_PART_PROP,partID) = real(blockID)
        particles(PROC_PART_PROP,partID) = real(pt_meshMe)
        particles(POSX_PART_PROP:POSZ_PART_PROP,partID) = pos(IAXIS:KAXIS)


        ! Now do velocities.  Originally this code (and FLASH2) had xvel=yvel=zvel=0
        ! Alan suggested that it should mimic what Lattice did, and initialize from the mesh
        call Particles_mapFromMesh(pt_typeInfo(PART_MAPMETHOD,myType),pt_velNumAttrib,pt_velAttrib,pos,&
             bndbox,delta,solnData,partAttribVec)
           
        do ind=1,pt_velNumAttrib
           particles(pt_velAttrib(PART_DS_IND,ind),partID) = partAttribVec(ind)
        end do
     endif
     
  enddo
  
  !---------------------------------------------------------------------
  ! Cleanup

  ! Release the mesh data pointer for this block.
  call Grid_releaseBlkPtr(blockID,solnData)

  ! Set the particle database local number of particles.

  pt_numLocal = pt_numLocal + pAdd

  ! cleanup
  deallocate(icell2ijk)
  deallocate(mass)
  deallocate(cumuMass)
  deallocate(xLeftCoord)
  deallocate(yLeftCoord)
  deallocate(zLeftCoord)

  return

end subroutine pt_initPositionsWithDensity


