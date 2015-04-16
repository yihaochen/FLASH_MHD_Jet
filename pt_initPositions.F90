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
      pt_meshNumProcs
       

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

  integer       :: i
  integer, save :: pt_ind
  logical       :: IsInBlock
  real, dimension(3) :: blockSize, blockCenter, del, pos
  integer       :: blockType,mapType
  !logical       :: IsInBlock,IsInSphere
  integer       :: part_props=NPART_PROPS
  integer       :: nozzle=1, pt_customNumAttrib=5
  integer,dimension(PART_ATTR_DS_SIZE, 5) :: pt_customAttrib
!----------------------------------------------------------------------

  pt_customAttrib(PART_DS_IND,1)=DEN0_PART_PROP
  pt_customAttrib(GRID_DS_IND,1)=DENS_VAR
  pt_customAttrib(PART_DS_IND,2)=DENS_PART_PROP
  pt_customAttrib(GRID_DS_IND,2)=DENS_VAR
  pt_customAttrib(PART_DS_IND,3)=MAGX_PART_PROP
  pt_customAttrib(GRID_DS_IND,3)=MAGX_VAR
  pt_customAttrib(PART_DS_IND,4)=MAGY_PART_PROP
  pt_customAttrib(GRID_DS_IND,4)=MAGY_VAR
  pt_customAttrib(PART_DS_IND,5)=MAGZ_PART_PROP
  pt_customAttrib(GRID_DS_IND,5)=MAGZ_VAR

  !       Initialization now done in Particles_init.
  
  !        Particle slot number
  !print *, "pt_initPositions"
  pt_ind = pt_numLocal

  ! Location of block faces
  call Grid_getblkPhysicalSize(blockID, blockSize)
  call Grid_getBlkCenterCoords(blockID, blockCenter)
  !call Grid_getDeltas(blockID, del)

  do i = 1, sim_ptInitNum
     call pt_getRandomPos(pos)

     ! Check if particle is in this block 
     isInBlock = (maxval(abs(pos-blockCenter) - 0.5*blockSize) < 0.0)
     !write(*,'(4es11.3)') maxval(abs(pos-blockCenter)), blockSize
     ! If it is, keep it; otherwise discard it.
      if (IsInBlock) then
        pt_ind = pt_ind + 1
        particles(BLK_PART_PROP,pt_ind) = real(blockID)
        particles(PROC_PART_PROP,pt_ind) = real(pt_meshMe)
        particles(POSX_PART_PROP:POSZ_PART_PROP,pt_ind)  = pos
        particles(TAG_PART_PROP,pt_ind)   = pt_ind + pt_meshMe*pt_meshNumProcs
        particles(TADD_PART_PROP,pt_ind)  = dr_initialSimTime
        particles(TAU_PART_PROP,pt_ind)   = 0.0
        !particles(VELX_PART_PROP,pt_ind)  = 0.0
        !particles(VELY_PART_PROP,pt_ind)  = 0.0
        !print *, "pt_initPositions, particle pos are ", particles(POSX_PART_PROP:POSZ_PART_PROP,pt_ind) 
        !write(*,'(3i4, 3es11.3)') pt_meshMe, pt_ind,i, pos
      endif




  enddo
  !delTheta = 2.0*PI / sim_ptnumTheta

  !rxvec = cross(sim(nozzle)%jetvec, (/ 0.0, sim_smallx, 1.0 /))
  !rxvec = rxvec / sqrt(sum(rxvec(:)*rxvec(:)))
  !ryvec = cross(sim(nozzle)%jetvec, rxvec)
  !print *, "pt_initPositions begin loop and sim_ptnumTheta is", sim_ptnumTheta
  !do k = -1, 1, 2
  ! do j = 0, 1
  !    if (j == 0) then
  !       imax = 1
  !    else
  !       imax = sim_ptnumTheta
  !    endif

  !    do i = 0, imax-1
  !       pos = sim(nozzle)%pos+sim_smallx&
  !             + sim(nozzle)%jetvec*sim(nozzle)%length*0.9999*k&
  !             + sim(nozzle)%radius*0.9999*j*(rxvec*cos(delTheta*i) + ryvec*sin(delTheta*i))
  !       ! Check if particle is in this block 
  !       isInBlock = (maxval(abs(pos-blockCenter) - 0.5*blockSize) < 0.0)
  !       !write(*,'(4es11.3)') maxval(abs(pos-blockCenter)), blockSize
  !       ! If it is, keep it; otherwise discard it.
  !        if (IsInBlock) then
  !          pt_ind = pt_ind + 1
  !          particles(BLK_PART_PROP,pt_ind) = real(blockID)
  !          particles(PROC_PART_PROP,pt_ind) = real(pt_meshMe)
#ifdef MASS_PART_PROP
  !          particles(MASS_PART_PROP,pt_ind) = 1.
#endif
  !          particles(POSX_PART_PROP:POSZ_PART_PROP,pt_ind)  = pos
  !          particles(TAG_PART_PROP,pt_ind)   = real((k+1)/2*(sim_ptnumTheta+1)+j+i)
  !          !particles(VELX_PART_PROP,pt_ind)  = 0.0
  !          !particles(VELY_PART_PROP,pt_ind)  = 0.0
  !          !print *, "pt_initPositions, particle pos are ", particles(POSX_PART_PROP:POSZ_PART_PROP,pt_ind) 
  !          write(*,'(5i4, 3es11.3)') pt_meshMe, pt_ind,i,j,k, pos
  !        endif

  !    enddo
  ! enddo
  !enddo  


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
       pt_posAttrib,pt_customNumAttrib,pt_customAttrib,mapType)

  success=.true.


  return

!----------------------------------------------------------------------
  
end subroutine pt_initPositions


