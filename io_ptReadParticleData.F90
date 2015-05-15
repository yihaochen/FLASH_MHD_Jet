!!****if* source/IO/IOParticles/hdf5/serial/io_ptReadParticleData
!!
!! NAME
!!
!! io_ptReadParticleData
!!
!!
!! SYNOPSIS
!!
!! io_ptReadParticleData()
!!
!!
!!
!! DESCRIPTION
!!
!!    This routine reads out the particle data in a separate hdf5 file
!!    It calls  ioh5_read_particles
!!    This is being done to make particles easy to debug and since the particles
!!    are not associated with the mesh data it makes since to separate it from
!!    the checkpoint files and plotfiles
!!
!! ARGUMENTS
!!
!!
!!
!! NOTES
!!
!!
!!***


subroutine io_ptReadParticleData()

  use IO_data, ONLY : io_outputSplitNum, &
       io_chkptFileID, io_globalMe, io_globalNumProcs
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_getPrev, &
    RuntimeParameters_get
  use Particles_interface, ONLY : Particles_putLocalNum
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_putLocalNumBlks

  use Grid_data, ONLY : gr_globalNumBlocks

  use Particles_data, ONLY : particles, pt_maxPerProc, pt_posInitialized
  use IO_interface, ONLY : IO_getScalar



  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"


  integer :: localNumParticles, ierr, i, particleOffset
  integer, save :: globalNumParticles !done for IBM compilers
  integer :: blocksPerFile, globalNumBlocks, blkOffset, localNumBlocks


  integer :: localNumBlockst, reLocalNumParticles, lb, startIndex, endIndex
  logical :: useParticles

  integer :: status(MPI_STATUS_SIZE), jproc
  
  real, allocatable :: particlest(:,:)

  integer :: particlesPerBlk(MAXBLOCKS)

    !for property-by-property read-in.
  integer :: fileNumPartProps, j, maxPartProps
  character(len=24), allocatable :: filePropNames(:) 
!!$  integer :: fileToCurrentMap(NPART_PROPS)
  integer,allocatable :: fileToCurrentMap(:)
  integer :: propIndex
  character(len=24) :: propString

  !It is necessary to initialize localNumParticles to 0 for the 
  !Particles_putLocalNum subroutine call because localNumParticles
  !will only get a new value when an MPI rank owns at least 1 block.
  localNumParticles = 0

  pt_posInitialized = .false. !So it is set even if we return early.

  !allocate particles data structure
  !need to get previous runtime parameter storing max particles per proc
  call RuntimeParameters_getPrev("pt_maxPerProc", pt_maxPerProc)

  !see if we are even running particles, if not, then return
  call RuntimeParameters_get("useParticles", useParticles)

  if(.not. useParticles) then
     return
  end if

  !Added call to RuntimeParameter to ensure that proc 0 gets the correct maximum
  !number of particles per block, if on a restart the number changes.  Appears
  !only to be an issue on gin and mongchi.
  call RuntimeParameters_get("pt_maxPerProc", pt_maxPerProc)

  allocate (particles(NPART_PROPS,pt_maxPerProc), stat=ierr)
  if (ierr /= 0) then
     call Driver_abortFlash("io_readParticlesData:  could not allocate particle array")
  endif

  !particles must be initialized or the entire particles algorithm will fail
  particles = NONEXISTENT

  !get the global number of particles in the simulation
  call IO_getScalar("globalNumParticles", globalNumParticles)

  if(globalNumParticles <= 0 ) then
     return
  else
     pt_posInitialized = .true.
  end if

  !allocate space for the temporary particles data struct
  allocate (particlest(NPART_PROPS,pt_maxPerProc), stat=ierr)
  if (ierr /= 0) then
     call Driver_abortFlash("io_readParticlesData:  could not allocate particlet array")
  endif


  
  !in UG this is always 1
  call Grid_getLocalNumBlks(localNumBlocks)
     

  particleOffset = 0
  blkOffset = 0

  call Grid_putLocalNumBlks(localNumBlocks)
  

  !in this serial implementation only the master processor will read all the data
  if (io_globalME == MASTER_PE) then
     
     do jproc=0, io_globalNumProcs -1
   
        if (jproc /= MASTER_PE) then
           
           call MPI_RECV(localNumBlockst, 1, FLASH_INTEGER, & 
                jproc, 1, MPI_COMM_WORLD, status, ierr)
           
        else
           
           localNumBlockst = localNumBlocks

        end if

        if (localNumBlockst > 0) then

           !return an array particlePerBlk holding the number 
           !of particles on each blk on the local proc
           call io_h5read_localnp(io_chkptFileID, &
                particlesPerBlk, &
                localNumBlockst, &
                gr_globalNumBlocks, &
                blkOffset)
        
                      
           !now find the newLocalNumParticles
           reLocalNumParticles = 0
           do lb=1, localNumBlockst
              reLocalNumParticles = particlesPerBlk(lb) + reLocalNumParticles
           end do
           
           
           if (reLocalNumParticles > pt_maxPerProc) then
              call Driver_abortFlash &
                   ('[io_ptReadParticleData] ERROR: too many particles on this proc; increase pt_maxPerProc')
           end if
           
           call io_h5read_num_props(io_chkptFileID, fileNumPartProps)
           !print *, "NumPartProps: ",fileNumPartProps
           
           allocate(filePropNames(fileNumPartProps))
           maxPartProps = MAX(fileNumPartProps, NPART_PROPS)
           allocate(fileToCurrentMap(maxPartProps))
           fileToCurrentMap = NONEXISTENT
           
           
           
           call io_h5read_particle_names(io_chkptFileID, filePropNames, fileNumPartProps);
           
           
          ! print *, filePropNames
           !generate mapping
           
           do i = 1,NPART_PROPS
              
              !iterate over part
              call Simulation_mapIntToStr(i, propString, MAPBLOCK_PART)
              do j = 1, fileNumPartProps
                 
                 if(propString .eq. filePropNames(j)) then
                   ! print *, propString, filePropNames(j)
                    fileToCurrentMap(j) = i
                    exit
                 end if
                 !print *, propString, filePropNames(i)
                 
              end do
              
              
           end do
           
           deallocate(filePropNames)
           !call Driver_abortFlash("DEBUG ABORT")
           
           do i = 1, fileNumPartProps
              
              if(fileToCurrentMap(i) .eq. NONEXISTENT ) then
                 !force iteration
                 cycle
              end if
              
           !   print *, "prop map:", i, fileToCurrentMap(i)
           !   print *, fileNumPartProps, NPART_PROPS
              call io_h5read_single_part_prop(io_chkptFileID, &
                   particlest, &
                   reLocalNumParticles, &
                   fileToCurrentMap(i), &
                   i, &
                   fileNumPartProps, &
                   NPART_PROPS, &
                   particleOffset)

            !  print *, particles(i,:)
           end do
           deallocate(fileToCurrentMap)
           
!DEV: Hold onto this for time being if reads need to be sped up --PR         
!!$           call io_h5read_particles(io_chkptFileID, &
!!$                particlest, &
!!$                reLocalNumParticles, &
!!$                NPART_PROPS, &
!!$                particleOffset);           
           
           !reset particles BLK_PART_PROP because it could have changed on restart
           startIndex = 1
           do lb=1, localNumBlockst
              
              if(particlesPerBlk(lb) > 0) then
                 endIndex = startIndex + particlesPerBlk(lb)
                 
                 particlest(BLK_PART_PROP,startIndex:endIndex-1) = lb
                 startIndex = endIndex
                 
              end if
           end do
           
        
!!$        !read the local number of particles on each blk/proc
!!$        !only relevant in UG because 1 blk per proc.  This data struct
!!$        !only exists in UG
!!$        call io_h5read_localnp(io_chkptFileID, &
!!$             localNumParticlest, &
!!$             localNumBlocks, &
!!$             gr_globalNumBlocks, &
!!$             blkOffset)
!!$
!!$        !read particles into temp particles datastruct
!!$        call io_h5read_particles(io_chkptFileID, &
!!$             particlest, &
!!$             localNumParticlest, &
!!$             NPART_PROPS, &
!!$             particleOffset)

           !increment the particle Offset
           particleOffset = particleOffset + relocalNumParticles
           blkOffset = blkOffset + localNumBlockst


           !now the MASTER_PE sends the data to the other procs
           if (jproc /= MASTER_PE) then

              call MPI_SEND (relocalNumParticles, &
                   1, & 
                   FLASH_INTEGER, jproc, 2, MPI_COMM_WORLD, ierr)
              
              
              call MPI_SEND (particlest(1,1), &
                   relocalNumParticles*NPART_PROPS, & 
                   FLASH_REAL, jproc, 3, MPI_COMM_WORLD, ierr)
              
           else
              
              localNumParticles = relocalNumParticles
              particles(:,1:localNumParticles) = particlest(:,1:localNumParticles)
              
           end if !if jproc /= MASTER_PE
        end if  !if localNumBlockst > 0

     end do !jprocs = 0, io_globalNumProcs -1
     
  else !myPE /= MASTER_PE

     !now send localNumBlocks
     call MPI_SEND(localNumBlocks, 1, FLASH_INTEGER, & 
          MASTER_PE, 1, MPI_COMM_WORLD, ierr)

     if(localNumBlocks > 0) then

        !now other procs recv the data
        call MPI_RECV(localNumParticles, &
             1, & 
             FLASH_INTEGER, MASTER_PE, 2, & 
             MPI_COMM_WORLD, status, ierr)
        
        
        call MPI_RECV(particles(1,1), &
             localNumParticles*NPART_PROPS, & 
             FLASH_REAL, MASTER_PE, 3, & 
             MPI_COMM_WORLD, status, ierr)
        
     end if !localNumBlocks > 0
  end if !(myPE == MASTER_PE)
  

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  call Particles_putLocalNum(localNumParticles)
  deallocate(particlest)

  return

end subroutine io_ptReadParticleData
