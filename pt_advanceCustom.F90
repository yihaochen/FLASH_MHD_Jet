!!****if* source/Particles/localAPI/pt_advanceCustom
!!
!! NAME
!!
!!  pt_advanceCharged
!!
!! SYNOPSIS
!!
!!  call pt_advanceCharged(real(in)   :: dtOld,
!!                        real(in)    :: dtNew,
!!                        real(inout) :: particles(:,p_count),
!!                        integer(in) :: p_count)
!!                        integer(in) :: ind)
!!
!! DESCRIPTION
!!
!!   Advances particles in time, based on rungekutta (pt_advanceRK.F90).
!!
!! ARGUMENTS
!!
!!   dtOld : previous time interval 
!!   dtNew : current time interval
!!   particle -- particles on which to operate
!!   p_count - the number of particles in the list to advance
!!   ind -- index into pt_typeInfo
!!
!! NOTES
!!
!!  No special handling is done for the first call - it is assumed that particle
!!  initialization fills in initial velocity components properly.
!!
!!
!!***

subroutine pt_advanceCustom(dtOld,dtNew, particles,p_count, ind)
    
  
  use Particles_data, ONLY: pt_numLocal, pt_maxPerProc, &
       useParticles, pt_typeInfo, &
       pt_gcMaskForAdvance, pt_gcMaskSizeForAdvance, pt_meshMe, &
       pt_posAttrib, pt_velNumAttrib,pt_velAttrib, pt_meshNumProcs

  use Driver_data, ONLY : dr_simTime, dr_initialSimTime
  use Hydro_data, ONLY : hy_bref
  use Grid_interface, ONLY : Grid_mapMeshToParticles
  use Simulation_data, ONLY : sim_ptAddPeriod
  
  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Particles.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  integer, INTENT(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(INOUT) :: particles

  integer       :: i,particleTypes
  
  integer,dimension(MAXBLOCKS, 1) :: perBlk
  real          :: jumpx,jumpy,jumpz
  real,allocatable :: origVel(:,:)
  integer :: part_props=NPART_PROPS
  integer :: mapType 
  
  integer      :: pt_customNumAttrib=4
  integer,dimension(PART_ATTR_DS_SIZE, 4) :: pt_customAttrib
  integer      :: pt_newParticleNumAttrib=5
  integer,dimension(PART_ATTR_DS_SIZE, 5) :: pt_newParticleAttrib
  integer              :: clock
  integer,dimension(2) :: seed
  real  :: prob, rho13, A
  real, dimension(MDIM,1) ::  pos
  logical :: addNewSuccess
!!------------------------------------------------------------------------------
  
  ! Don't do anything if runtime parameter isn't set
  if (.not.useParticles ) return


  ! Some initializations, might need to move them to Particle_init
  pt_newParticleAttrib(PART_DS_IND,1)=DEN0_PART_PROP
  pt_newParticleAttrib(GRID_DS_IND,1)=DENS_VAR
  pt_newParticleAttrib(PART_DS_IND,2)=DENS_PART_PROP
  pt_newParticleAttrib(GRID_DS_IND,2)=DENS_VAR
  pt_newParticleAttrib(PART_DS_IND,3)=MAGX_PART_PROP
  pt_newParticleAttrib(GRID_DS_IND,3)=MAGX_VAR
  pt_newParticleAttrib(PART_DS_IND,4)=MAGY_PART_PROP
  pt_newParticleAttrib(GRID_DS_IND,4)=MAGY_VAR
  pt_newParticleAttrib(PART_DS_IND,5)=MAGZ_PART_PROP
  pt_newParticleAttrib(GRID_DS_IND,5)=MAGZ_VAR

  pt_customAttrib(PART_DS_IND,1)=DENS_PART_PROP
  pt_customAttrib(GRID_DS_IND,1)=DENS_VAR
  pt_customAttrib(PART_DS_IND,2)=MAGX_PART_PROP
  pt_customAttrib(GRID_DS_IND,2)=MAGX_VAR
  pt_customAttrib(PART_DS_IND,3)=MAGY_PART_PROP
  pt_customAttrib(GRID_DS_IND,3)=MAGY_VAR
  pt_customAttrib(PART_DS_IND,4)=MAGZ_PART_PROP
  pt_customAttrib(GRID_DS_IND,4)=MAGZ_VAR

  mapType=pt_typeInfo(PART_MAPMETHOD,ind)

!!------------------------------------------------------------------------------
  ! Update the particle positions to temporary ("predicted") values
  do i = 1, p_count
 
     jumpx = dtNew * particles(VELX_PART_PROP,i)
     particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + jumpx

     if(NDIM >1) then
        jumpy = dtNew * particles(VELY_PART_PROP,i)
        particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + jumpy
        if(NDIM >2) then
           jumpz = dtNew * particles(VELZ_PART_PROP,i)
           particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + jumpz
        end if
     end if

  enddo


  ! Now save the original velocity values
  allocate(origVel(p_count,MDIM))
  origVel(:,1) = particles(VELX_PART_PROP,1:p_count)
  origVel(:,2) = particles(VELY_PART_PROP,1:p_count)
  origVel(:,3) = particles(VELZ_PART_PROP,1:p_count)

  ! Map the updated gas velocity field at the temporary positions to
  ! obtain a second estimate of velocities;

  call Grid_mapMeshToParticles(particles,&
       part_props, BLK_PART_PROP, p_count,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)

  ! Adjust particle positions, using the second point velocities
  do i = 1, p_count
     particles(POSX_PART_PROP,i) =  particles(POSX_PART_PROP,i) + &
          dtNew * 0.5*(particles(VELX_PART_PROP,i) - origVel(i,1))
     if(NDIM>1)&
          particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + &
          dtNew * 0.5*(particles(VELY_PART_PROP,i) - origVel(i,2))
     if(NDIM>2)&
          particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + &
          dtNew * 0.5*(particles(VELZ_PART_PROP,i) - origVel(i,3))
  enddo

  ! done with this temporary storage
  deallocate(origVel)



  ! Map the updated gas velocity field onto the current particle positions to
  ! obtain the updated particle velocities - for the next integration step
  ! as well as for particle plot files etc.

  call Grid_mapMeshToParticles(particles,&
       part_props, BLK_PART_PROP,p_count,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)

  ! Map the custom fields onto the partitcles
  call Grid_mapMeshToParticles(particles,&
       part_props,BLK_PART_PROP, p_count,&
       pt_posAttrib,pt_customNumAttrib,pt_customAttrib,mapType)

  ! update the synchrotron lifetime and cutoff gamma
  do i = 1, p_count
    rho13 = (particles(DENS_PART_PROP,i)/particles(DEN0_PART_PROP,i))**(1.0/3.0)
    ! e^4/(me^3*c^5) = 2.907728E-9
    ! if units == "none" (default) -> hy_bref=1.0
    ! if units == "cgs" -> hy_bref=sqrt(4*pi)
    A = 4.0/9.0*2.907728E-9*4.0*PI/hy_bref/hy_bref&
        *sum(particles(MAGX_PART_PROP:MAGZ_PART_PROP,i)*particles(MAGX_PART_PROP:MAGZ_PART_PROP,i))
        
    particles(TAU_PART_PROP,i) = particles(TAU_PART_PROP,i) + rho13*A*dtNew
    particles(GAMC_PART_PROP,i) = rho13 / particles(TAU_PART_PROP,i)

  enddo

!!------------------------------------------------------------------------------
  ! Add new particles stochastically

  !write(*,'(i5, A28, i5)') pt_meshMe, 'Before addNew, pt_numLocal=', pt_numLocal
  !write(*,'(i5, A28, i5)') pt_meshMe, 'Before addNew, p_count    =', p_count
  call RANDOM_NUMBER(prob)
  !write(*,'(i5, 2f9.5)') pt_meshMe, prob, 1.0/pt_meshNumProcs/sim_ptAddPeriod*dtNew
  !write(*,'(i5, f9.5, es11.3)') pt_meshMe, prob, dtNew
  if (prob .le. 1.0/pt_meshNumProcs/sim_ptAddPeriod*dtNew) then
     call pt_getRandomPos(pos(:,1))
     !write(*,'(i5, A25, 3es11.3)') pt_meshMe, 'Adding a new particle at', pos(:,1)
     call Particles_addNew(1, pos, addNewSuccess) 
     if (addNewsuccess) then
        !write(*,'(i5, A25, 3es11.3)') pt_meshMe, 'Added a new particle at', pos(:,1)
     end if
  else
      call Particles_addNew(0, pos, addNewSuccess)
  end if

  !write(*,'(i5, A28, i5)') pt_meshMe, 'After addNew, pt_numLocal=', pt_numLocal
  !write(*,'(i5, A28, i5)') pt_meshMe, 'After addNew, p_count    =', p_count
  do i = 1, pt_numLocal
  ! if the particle is newly added, map the properties first
     if (particles(DEN0_PART_PROP,i) .le. 0.0 ) then
        call Grid_mapMeshToParticles(particles(:,i:i),&
             part_props,BLK_PART_PROP, 1,&
             pt_posAttrib,pt_newParticleNumAttrib,pt_newParticleAttrib,mapType)

        A = 4.0/9.0*1.938486E-09*4.0*PI/hy_bref**2&
            *sum(particles(MAGX_PART_PROP:MAGZ_PART_PROP,i)*particles(MAGX_PART_PROP:MAGZ_PART_PROP,i))
            
        particles(TAU_PART_PROP,i) = particles(TAU_PART_PROP,i) + A*dtNew
        particles(GAMC_PART_PROP,i) = 1.0 / particles(TAU_PART_PROP,i)
     end if
  enddo
  
end subroutine pt_advanceCustom
