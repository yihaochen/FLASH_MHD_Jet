!!****if* source/Simulation/SimulationMain/magnetoHD/MHD_Jet/
!!
!! NAME
!!    pt_initPositions
!!
!! SYNOPSIS
!!
!!    call pt_initPositions(integer(in)  :: blockID,
!!                          logical(OUT) :: success)
!!
!! DESCRIPTION
!!    Initialize particle locations.  This version sets up particles
!!      which are evenly distributed in circular space
!!
!! ARGUMENTS
!!
!!  blockID:        local block ID containing particles to create
!!  success:        returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!! PARAMETERS
!!
!!    pt_numX:      number of particles along physical x-axis of domain
!!    pt_numY:      number of particles along physical y-axis of domain
!!    pt_numZ:      number of particles along physical z-axis of domain
!!
!!***


subroutine pt_initPositions (blockID,success)

!  use Particles_data

  use Particles_data, ONLY:  pt_numLocal, particles,&
      pt_posAttrib,pt_velNumAttrib, pt_velAttrib,pt_typeInfo, pt_meshMe,&
      pt_meshNumProcs, pt_newParticleNumAttrib, pt_newParticleAttrib,&
      pt_numParticlesWanted


  use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_mapMeshToParticles, &
                             Grid_getBlkPhysicalSize, Grid_getBlkCenterCoords, &
                             Grid_getDeltas
  use Driver_data, ONLY : dr_initialSimTime

  use Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  integer, INTENT(in) :: blockID
  logical,intent(OUT) :: success

  integer       :: i, k, nozzle=1
  integer, save :: pt_ind
  logical       :: IsInBlock
  real, dimension(3) :: blockSize, blockCenter, del, pos
  integer       :: blockType,mapType
  !logical       :: IsInBlock,IsInSphere
  integer       :: part_props=NPART_PROPS
!----------------------------------------------------------------------

  !       Initialization now done in Particles_init.


  !print *, "pt_initPositions"
  !        Particle slot number
  pt_ind = pt_numLocal

  ! Location of block faces
  call Grid_getblkPhysicalSize(blockID, blockSize)
  call Grid_getBlkCenterCoords(blockID, blockCenter)
  !call Grid_getDeltas(blockID, del)

  !do i = 1, sim_ptInitNum
  !   call pt_getRandomPos(pos)

  !   ! Check if particle is in this block
  !   isInBlock = (maxval(abs(pos-blockCenter) - 0.5*blockSize) < 0.0)
  !   !write(*,'(4es11.3)') maxval(abs(pos-blockCenter)), blockSize
  !   ! If it is, keep it; otherwise discard it.
  !   if (IsInBlock) then
  !      pt_ind = pt_ind + 1
  !      particles(BLK_PART_PROP,pt_ind) = real(blockID)
  !      particles(PROC_PART_PROP,pt_ind) = real(pt_meshMe)
  !      particles(POSX_PART_PROP:POSZ_PART_PROP,pt_ind)  = pos
  !      particles(TAG_PART_PROP,pt_ind)   = pt_ind + pt_meshMe*pt_meshNumProcs
  !      particles(TADD_PART_PROP,pt_ind)  = dr_initialSimTime
  !      particles(TAU_PART_PROP,pt_ind)   = 0.0
  !      !particles(VELX_PART_PROP,pt_ind)  = 0.0
  !      !particles(VELY_PART_PROP,pt_ind)  = 0.0
  !      !print *, "pt_initPositions, particle pos are ", particles(POSX_PART_PROP:POSZ_PART_PROP,pt_ind)
  !      !write(*,'(3i4, 3es11.3)') pt_meshMe, pt_ind,i, pos
  !   endif
  !enddo

  !print *, "pt_initPositions begin loop"
  do k = -1, 1, 2
     pos = sim(nozzle)%pos+sim_smallx*sim(nozzle)%radius&
           + sim(nozzle)%jetvec*sim(nozzle)%length*0.9999*k
     ! Check if particle is in this block
     isInBlock = (maxval(abs(pos-blockCenter) - 0.5*blockSize) < 0.0)
     !write(*,'(i4, es11.3)') pt_meshMe, maxval(abs(pos-blockCenter) - 0.5*blockSize)
     ! If it is, keep it; otherwise discard it.
     if (IsInBlock) then
       pt_ind = pt_ind + 1
       particles(BLK_PART_PROP,pt_ind) = real(blockID)
       particles(PROC_PART_PROP,pt_ind) = real(pt_meshMe)
       particles(POSX_PART_PROP:POSZ_PART_PROP,pt_ind)  = pos
       if (k.eq.-1) then
          particles(TAG_PART_PROP,pt_ind) = pt_numParticlesWanted + 1
       else
          particles(TAG_PART_PROP,pt_ind) = pt_numParticlesWanted + 2
       endif
       particles(TADD_PART_PROP,pt_ind)  = dr_initialSimTime
       particles(TAU_PART_PROP,pt_ind)   = 0.0
       !particles(VELX_PART_PROP,pt_ind)  = 0.0
       !particles(VELY_PART_PROP,pt_ind)  = 0.0
       !print *, "pt_initPositions, particle pos are ", particles(POSX_PART_PROP:POSZ_PART_PROP,pt_ind)
     endif

  enddo


  !       Setting the particle database local number of particles
  pt_numLocal = pt_ind
  !print *, "pt_numLocal", pt_numLocal
  !       Now initialize velocity properties for the new particles
  mapType=pt_typeInfo(PART_MAPMETHOD,1)
  call Grid_mapMeshToParticles(particles,&
       part_props,BLK_PART_PROP, pt_numLocal,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)

  ! Map the custom fields onto the partitcles
  call Grid_mapMeshToParticles(particles,&
       part_props,BLK_PART_PROP, pt_numLocal,&
       pt_posAttrib,pt_newParticleNumAttrib,pt_newParticleAttrib,mapType)

  success=.true.


  return

!----------------------------------------------------------------------

end subroutine pt_initPositions


