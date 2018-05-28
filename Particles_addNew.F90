!!****if* source/Particles/ParticlesInitialization/Particles_addNew
!!
!! NAME
!!    Particles_addNew
!!
!! SYNOPSIS
!!    call Particles_addNew( integer(in)  :: count,
!!                  optional,real(in)     :: pos(MDIM,count),
!!                           logical(out) :: success)
!!
!! DESCRIPTION
!!
!!    This routine allows particles to be added during evolution.
!!    In the particles data structure it always initializes the tag and
!!    processor ID. If the optional argument "pos" is present then it
!!    will also initialize the position and block ID attributes in the
!!    particles. It returns the value FALSE in success if there isn't
!!    enough space left in the data structure for the requested number
!!    of particles.
!!
!! ARGUMENTS
!!
!!     count   :: the count of particles to be added
!!     pos     :: optional, contains the coordinates of the particles
!!     success :: This arg returns TRUE if there was enough space 
!!                in the particles data structure to add the requested
!!                number of particles, FALSE otherwise.
!!
!!    
!!  NOTES
!!
!!   This routine must be called collectively, i.e., by all MPI tasks
!!   in the pt_meshComm communicator, if the optional "pos" argument
!!   is present.
!!
!!   The constant MDIM is defined in constants.h .
!!
!!***

!!#define DEBUG_PARTICLES

subroutine Particles_addNew (count, pos, shock, success)
  
  use Particles_data, ONLY : particles, &
       pt_maxPerProc, pt_numLocal, pt_meshComm, pt_meshMe, pt_indexList, &
       pt_indexCount, useParticles

  use Logfile_interface, ONLY : Logfile_stampMessage
  use Grid_interface, ONLY : Grid_moveParticles
  use pt_interface, ONLY : pt_findTagOffset
  use Driver_data, ONLY : dr_simTime
  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Particles.h"
#include "Flash_mpi.h"

  integer, INTENT(in) :: count
  real, optional, dimension(MDIM,count), intent(IN)::pos
  real, optional, intent(IN):: shock
  logical, intent(OUT) :: success

  integer :: i, tagOffset, ierr, effCount
  logical,parameter :: coords_in_blk=.true.
  logical :: doAll, doAny, doLocal(2),doGlobal(2)
  character(len=80) :: message

  !integer      :: mapType=PART_MAPMETHOD

!-------------------------------------------------------------------------------
  if (.not.useParticles) return

  doLocal(1) = ((pt_numLocal+count).le.pt_maxPerProc)
  doLocal(2) = (.NOT. doLocal(1))

  !write(*,*) pt_meshMe, 'AllReduce'
  call MPI_AllReduce(doLocal(1),doGlobal(1),2,MPI_LOGICAL,MPI_LAND,pt_meshComm,ierr)
  doAll = doGlobal(1)
  doAny = .NOT. doGlobal(2)

!  if((pt_numLocal+count).le.pt_maxPerProc) then
  if(doAny) then
     !write(*,*) pt_meshMe, 'addNew - doAny'
     if (doLocal(1)) then
        !write(*,*) pt_meshMe, 'addNew - doAny - doLocal', count
        success=.true.
        effCount = count
     else
        !write(*,*) pt_meshMe, 'addNew - doAny - no local', count
        success = (count==0)
        effCount = 0
     end if
     !write(*,*) pt_meshMe, 'findTagOffset, effCount=', effCount
     call pt_findTagOffset(effCount,tagOffset)
     !write(*,*) pt_meshMe, 'findTagOffset, tagOffset=', tagOffset
     do i = 1,effCount
        particles(PROC_PART_PROP,pt_numLocal+i) = pt_meshMe
        particles(TAG_PART_PROP,pt_numLocal+i)  = tagOffset+i
        particles(TADD_PART_PROP,pt_numLocal+i) = dr_simTime
        ! Set den0 < 0.0 for initialization in pt_advanceCustom
        particles(DEN0_PART_PROP,pt_numLocal+i) = -1.0
        particles(DENS_PART_PROP,pt_numLocal+i) = -1.0
        call pt_resetShockVars(particles(:,pt_numLocal+i), 1, 100.0, dr_simTime)
        call pt_resetShockVars(particles(:,pt_numLocal+i), 2, 100.0, dr_simTime)
        call pt_resetShockVars(particles(:,pt_numLocal+i), 3, 100.0, dr_simTime)
        particles(TAU0_PART_PROP,pt_numLocal+i) = 1E-100
        particles(CMB0_PART_PROP,pt_numLocal+i) = 1E-100
        particles(ICT0_PART_PROP,pt_numLocal+i) = 1E-100
        ! Which shock is this particle located
        ! Set to a large initial value, i.e. very weak shock
        particles(WHCH_PART_PROP,pt_numLocal+i)  = 100.0
        particles(GAMC_PART_PROP,pt_numLocal+i)  = 1E100
        particles(TYPE_PART_PROP,pt_numLocal+i) = PASSIVE_PART_TYPE
        if (present(shock)) then
           particles(SHOK_PART_PROP,pt_numLocal+i) = shock
        endif
     end do
     if(present(pos)) then
        !write(*,*) pt_meshMe, 'present(pos)'
        do i = 1,effCount
           particles(POSX_PART_PROP:POSZ_PART_PROP,pt_numLocal+i)=&
                pos(IAXIS:KAXIS,i)
        end do
        
        pt_numLocal=pt_numLocal+effCount
        !write(*,*) pt_meshMe, 'moveParticles'
        call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,&
             pt_numLocal,pt_indexList,pt_indexCount,coords_in_blk)

        !write(*,*) pt_meshMe, 'moveParticles - finished'
     else
        pt_numLocal=pt_numLocal+effCount
        !write(*,*) pt_meshMe, 'else moveParticles'
        call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,&
             pt_numLocal,pt_indexList,pt_indexCount,coords_in_blk)

        !write(*,*) pt_meshMe, 'else moveParticles - finished'
     end if
     if (.NOT. success) then
98      format('Particles_addNew for',I6,' particles failed')
        write(message,98) count
        print *,trim(message),' on',pt_meshMe
        call Logfile_stampMessage(message, force=.TRUE.)
     end if
  else
99   format('Particles_addNew for',I6,' particles failed everywhere') 
     write(message,99) count
     print *,trim(message),', detected on',pt_meshMe
     call Logfile_stampMessage(message, force=.TRUE.)
     success=.false.
  end if
  
  return
  
end subroutine Particles_addNew
