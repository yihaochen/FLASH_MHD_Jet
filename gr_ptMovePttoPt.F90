!!****if* source/Grid/GridParticles/GridParticlesMove/gr_ptMovePttoPt
!!
!! NAME
!!
!!  gr_ptMovePttoPt
!!
!! SYNOPSIS
!!
!!  gr_ptMovePttoPt (real    (inout) :: dataBuf (:,:),
!!                   integer (inout) :: localCount,
!!                   integer (in)    :: propCount,
!!                   integer (in)    :: maxCount,
!!                   integer (in)    :: numDest,
!!         optional, integet (in)    :: matchProp)
!!  
!! DESCRIPTION 
!!  
!!  This routine moves non stationary Grid data elements to the correct block and
!!  processor. The data elements such as particles may have moved because of time
!!  advancement, or after the grid has refined or derefined. Other data elements
!!  such as those associated with rays follow the path of the ray.
!!
!!  Overview of algorithm
!!
!!  * Sort the elements according to their target processors.
!!
!!  * Determine the number of elements to be sent to each processor
!!    from the current processor and share that info with all other
!!    processors.
!!
!!  * Post first all necessary non-blocking receives from the receiving
!!    processors using appropriate non-gapped memory locations for the
!!    receiving elements.
!!
!!  * Post next all necessary blocking sends from the current processor.
!!
!!  * Check for all receives to finish and update the new element count.
!!
!! ARGUMENTS 
!!
!!  dataBuf:     List of data elements. This is a two dimensional real array,
!!               the first dimension representing each particle's properties
!!               and the second dimension is the counting index.
!!
!!  localCount:  While coming in, it contains the current number of elements
!!               mapped to this processor. After all the data structure movement,
!!               the number of local elements might have change, and the new
!!               value is put back into it.
!!
!!  propCount:   Number of element attributes.
!!
!!  maxCount:    This is a parameter determined at runtime. It is the maximum
!!               number of elements that a simulation expects to hold in the
!!               domain at any time during a simulation. All the arrays in the
!!               elements unit are allocated based on this number.
!!
!!  numDest:     The count of data elements to be moved.
!!
!!  matchProp:   (Obsolete?) Property to be matched to find the destination of
!!               the data element
!!
!! NOTES
!!
!!  This is a cleanup and new rewrite of the original gr_ptMovePttoPt routine to
!!  avoid use of redundant memory during element movement and consequently avoiding
!!  code abortion due to larger memory request than number of elements in the
!!  domain at any time of the simulation __N.F. (May 2016).
!!
!! SEE ALSO
!!
!!  gr_ptLocalMatch
!!
!!***

subroutine gr_ptMovePttoPt (dataBuf, propCount, maxCount, localCount, numDest)

  use Timers_interface, ONLY : Timers_start, &
                               Timers_stop

  use Grid_data,        ONLY : gr_useParticles,       &
                               gr_meshNumProcs,       &
                               gr_meshMe,             &
                               gr_meshComm,           &
                               gr_useEnergyDeposition

  use gr_ptData,        ONLY : gr_ptBlkList,   &
                               gr_ptBlkCount,  &
                               gr_ptDestBuf,   &
                               gr_ptSourceBuf, &
                               gr_ptBlk,       &
                               gr_ptProc,      &
                               gr_ptNumToReduce

  use ut_sortInterface, ONLY : ut_sortOnProcs

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent (inout) :: localCount
  integer, intent (in)    :: propCount, maxCount, numDest
  real,    intent (inout) :: dataBuf (1:propCount, 1:maxCount)

  logical :: sendingData, receivingData

  integer :: sendCount, recvCount
  integer :: ierr
  integer :: numRecv
  integer :: proc
  integer :: sendParticles, recvParticles, procParticles
  integer :: gr_ptNumToReduce_old

  integer, parameter :: tagData = 33

  integer :: perProc   (1:gr_meshNumProcs)
  integer :: toProcs   (1:gr_meshNumProcs)
  integer :: fromProcs (1:gr_meshNumProcs)

  integer, allocatable :: req    (:)
  integer, allocatable :: status (:,:)
!
!
!     ...Immediate return, if no movement of elements needed.
!
!
  if (.not. (gr_useParticles .or. gr_useEnergyDeposition)) then
      return
  end if

  if (gr_meshNumProcs == 1) then
      return
  end if
!
!
!     ...Sort the particles to be moved based upon the processor number of 
!        their destination. Previously in the calling routine (Grid_moveParticles)
!        the dataBuf array has been copied to the gr_ptDestBuf array. It is
!        only passed in argument here to store the final new received particles
!        sitting in the gr_ptSourceBuf array.
!
!
  call ut_sortOnProcs (numDest,         &  ! number of particles currently in gr_ptDestBuf
                       propCount,       &  ! number of particle properties in gr_ptDestBuf
                       gr_ptProc,       &  ! processor info property index in gr_ptDestBuf
                       gr_meshNumProcs, &  ! total number of processors under consideration
                       gr_ptDestBuf,    &  ! particle array to be sorted
                       gr_ptSourceBuf,  &  ! temporary storage for particles while processing
                       perProc,         &  ! particle count destined for each processor
                       toProcs,         &  ! contains 1/0 (means send/no-send) to each processor
                       sendCount        )  ! number of processors that will send particles
!
!
!     ...Share that info with all other processors.
!
!
  call MPI_Alltoall (perProc,         & ! particle count destined for each processor on current processor
                     1,               & ! sending only 1 integer
                     FLASH_INTEGER,   & ! sending data type
                     fromProcs,       & ! particle count from each processor on current processor
                     1,               & ! receiving only 1 integer
                     FLASH_INTEGER,   & ! receiving data type
                     gr_meshComm,     & ! communicator
                     ierr             ) ! error handle
!
!
!     ...Check, if there is sufficient memory for consecutively storing the received
!        particles (without array gaps) on the current processor. Also check if the
!        dataBuf array will be able to hold all received particles.
!
!
  sendParticles = sum (perProc (:))
  recvParticles = sum (fromProcs (:))

  if (recvParticles > ubound (gr_ptSourceBuf,2)) then
      print *, "Overflow! Mesh PE", gr_meshMe, &
               ", receive particles count", recvParticles, &
               ", max particles for gr_ptSourceBuf", ubound (gr_ptSourceBuf,2)
      call Driver_abortFlash("[gr_ptMovePttoPt]: Insufficient space "//&
                             "in particles communication buffer: "//&
                             "increase pt_maxPerProc if you are using the Particles unit "//&
                             "or ed_maxRayCount if you are using the EnergyDeposition unit.")
  end if

  if (localCount + recvParticles > maxCount) then
      print *, "Overflow! Mesh PE", gr_meshMe, &
               ", receive particles total count", localCount + recvParticles, &
               ", max particle storage for dataBuf", maxCount
      ! Save the number
      gr_ptNumToReduce_old = gr_ptNumToReduce
      ! Remove at least the number of particles to allow receive from all procs
      gr_ptNumToReduce = max(localCount + recvParticles - maxCount, gr_ptNumToReduce)
      call gr_ptHandleExcess(dataBuf,propCount,localCount,maxCount)
      gr_ptNumToReduce = gr_ptNumToReduce_old

      !call Driver_abortFlash("[gr_ptMovePttoPt]: Insufficient space in particles data buffer.")
  end if
!
!
!     ...Post the non-blocking receives (if any) and the blocking sends (if any) to
!        move the particles between the processors.
!
!
  sendingData   = (sendParticles > 0)
  receivingData = (recvParticles > 0)

  if (receivingData) then
      numRecv = count (fromProcs (:) > 0)
      allocate (req (1:numRecv))
      numRecv = 0
      recvParticles = 0
      do proc = 0, gr_meshNumProcs - 1
         procParticles = fromProcs (proc + 1)
         if (procParticles > 0) then
             recvCount = propCount * procParticles
             call MPI_Irecv (gr_ptSourceBuf (1,recvParticles+1), &
                             recvCount,                          &
                             FLASH_REAL,                         &
                             proc,                               &
                             tagData,                            &
                             gr_meshComm,                        &
                             req (numRecv+1),                    &
                             ierr                                )

             numRecv = numRecv + 1
             recvParticles = recvParticles + procParticles
         end if
     end do
  end if

  if (sendingData) then
      sendParticles = 0
      do proc = 0, gr_meshNumProcs - 1
         procParticles = perProc (proc + 1)
         if (procParticles > 0) then
             sendCount = propCount * procParticles
             call MPI_Send (gr_ptDestBuf (1,sendParticles+1), &
                            sendCount,                        &
                            FLASH_REAL,                       &
                            proc,                             &
                            tagData,                          &
                            gr_meshComm,                      &
                            ierr                              )

             sendParticles = sendParticles + procParticles
         end if
     end do
  end if

  if (receivingData) then
      allocate (status (MPI_STATUS_SIZE, numRecv)) 
      call MPI_Waitall (numRecv, req, status, ierr)

      if (localCount + recvParticles > maxCount) then
         print *, "Overflow after receiving! Mesh PE", gr_meshMe, &
                  ", receive particles total count", localCount + recvParticles, &
                  ", max particle storage for dataBuf", maxCount
         ! Save the number
         gr_ptNumToReduce_old = gr_ptNumToReduce
         ! Remove at least the number of particles to allow receive from all procs
         gr_ptNumToReduce = max(localCount + recvParticles - maxCount, gr_ptNumToReduce)
         call gr_ptHandleExcess(dataBuf,propCount,localCount,maxCount)
         gr_ptNumToReduce = gr_ptNumToReduce_old
      end if

      dataBuf (1:propCount, localCount+1:localCount+recvParticles) = gr_ptSourceBuf (1:propCount, 1:recvParticles)
      localCount = localCount + recvParticles

      deallocate (req)
      deallocate (status)
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine gr_ptMovePttoPt
