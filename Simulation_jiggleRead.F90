!!****if* source/Simulation/SimulationMain/magnetoHD/MHD_Jet/Simulation_jiggleRead
!!
!! NAME
!!  Simulation_jiggleRead
!!
!! SYNOPSIS
!!
!!  use Simulation_jiggleRead
!!
!! DESCRIPTION
!!
!!  Reading the pointing director of the jet from a tabulated file and
!!  interpolate between time steps.
!!
!! ARGUMENTS
!!
!!
!!
!!
!!
!!***

subroutine Simulation_jiggleRead( nozzle, time, dt )

  use Simulation_data
  use Driver_data, ONLY : dr_globalMe, time, dt

  implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Simulation.h"

  integer, INTENT(in) :: nozzle
  real, INTENT(in) :: time, dt
  integer :: i, ierr, istat
  integer,dimension(1) :: i0arr
  integer,save :: i0
  character(len=250) :: line

  integer,save :: tableSize
  real,save,dimension(MAX_LINES_READ) :: table_time
  real,save,dimension(3,MAX_LINES_READ) :: table_jetvec
  real,save,dimension(3,MAX_LINES_READ) :: table_angvel
  real,dimension(3) :: jetaxis
  logical,save :: first_read=.true., eof_reached=.false.
  integer :: dummy

  if (first_read) then
     tableSize=MAX_LINES_READ
  else
     i0arr = minloc(abs(table_time(:)-time))
     i0 = i0arr(1)
  endif

  !if (dr_globalMe == MASTER_PE) then
  !   write(*,'(A3,i4,2e11.3)') 'i0:', i0, table_time(i0), table_time(i0+1)
  !endif

  do while ((i0.ge.tableSize).or.first_read)
  ! --------------------------------------------------------------------------
     ! Read the jet nozzle vectors from file.
     if (dr_globalMe == MASTER_PE) then
        i = 1
        open(sim_nozfileunit, file=sim_nozVecInput, status='OLD', iostat=istat)
        ! Skip the first line
        read(sim_nozfileunit,*)
        do while (i .le. tableSize)
           read(sim_nozfileunit, '(A)', iostat=istat) line

           ! Leave loop if end-of-file is reached
           if (istat < 0) then
              eof_reached = .true.
              tableSize=i-1
              write(*,*) 'eof_reached, tableSize', tableSize
              exit ! exit line reading
           endif

           ! Skip comment lines
           if (line(1:1) == '#') then
              cycle
           else
              read(line, 12) dummy, table_time(i), table_jetvec(:,i), table_angvel(:,i), dummy
12            format (1X, 1(I10, :, 1X), 7(es25.18, :, 1X), 1(I10, :, 1X))
              ! Start reading at the line of the file that corresponds to current
              ! simulation time
              ! TODO: need a better way to set the initial time
              if (table_time(i).lt.time-10.*dt) cycle
              i = i+1
           endif
        enddo
        close(sim_nozfileunit)
        write(*,'(A34, 2es11.3)') 'Read nozzle vectors between t = ', table_time(1), table_time(tableSize)
     endif

     ! Broadcast to other processors
     call MPI_BCAST(table_time,tableSize, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(table_jetvec,tableSize*3, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(table_angvel,tableSize*3, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

     first_read = .false.
     i0arr = minloc(abs(table_time(:)-time))
     i0 = i0arr(1)
     !if (dr_globalMe == MASTER_PE) then
     !   write(*,*) 'i0=', i0, 'tableSize=', tableSize
     !endif
     if (eof_reached.and.(i0.ge.tableSize)) then
        call Driver_abortFlash('[Simulation_jiggleRead] ERROR: simulation time out of nozzleVecInput scope')
     endif
  enddo
  ! --------------------------------------------------------------------------

  if ((table_time(i0)-time).gt.0.0) then
     i0 = i0-1
  endif

  !if (dr_globalMe == MASTER_PE) then
  !   write(*,*) 'i0:', i0, time, table_time(i0), table_time(i0+1)
  !   write(*,*) 'jetvec:', table_jetvec(:,i0), table_jetvec(:,i0+1)
  !endif

  jetaxis(:) = table_jetvec(:,i0) + &
  (table_jetvec(:,i0+1)-table_jetvec(:,i0))*(time-table_time(i0))/(table_time(i0+1)-table_time(i0))

  sim(nozzle)%angVel(:)=cross(sim(nozzle)%jetvec(:),jetaxis(:))/dt

  sim(nozzle)%jetvecOld(:)=sim(nozzle)%jetvec(:)
  sim(nozzle)%jetvec(:)=jetaxis(:)


end subroutine Simulation_jiggleRead
