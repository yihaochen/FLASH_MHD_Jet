!!****f* source/physics/sourceTerms/Heat/Heat
!!
!! NAME
!!
!!  Heat
!!
!!
!! SYNOPSIS
!!
!!  call Heat (integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount),
!!             real(IN)    :: dt,
!!             real(IN)    :: time)
!!
!!
!!
!! DESCRIPTION
!!
!!  Apply the stat+gauss source term operator to a block of
!!  zones. The energy generation rate is used to update the
!!  internal energy in the zone. The phonomenological heating
!!  rate is described as a 3-D Gauss function.
!!
!!  After we call stat+gauss, call the eos to update the
!!  pressure and temperature based on the phenomenological
!!  heating.
!!
!!
!! ARGUMENTS
!!
!!  blockCount : number of blocks to operate on
!!  blockList  : list of blocks to operate on
!!  dt         : current timestep
!!  time       : current time
!!
!!***

subroutine Heat (blockCount,blockList,dt,time)
!
!==============================================================================
!
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getDeltas, Grid_fillGuardCells
  use hy_uhd_interface, ONLY : hy_uhd_staggeredDivb
  !use Hydro_data, ONLY: hy_unsplitEosMode
  !use Eos_interface, ONLY : Eos_wrapped
  use Driver_data, ONLY : dr_globalMe, dr_nStep
  use Simulation_data
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Heat_data, ONLY : nPtProc, pos
  use Particles_data, ONLY : pt_randSeed
  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN)::blockList
  real,intent(IN) :: dt,time
  real, dimension(MDIM) :: del
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  logical :: gcMask(NUNK_VARS)

  integer :: blockID, blkInd, nozzle=1

  integer :: nPtNoz
  real,allocatable,dimension(:,:) ::  posNoz,pos_tmp
  real    :: prob
  logical      :: addNewSuccess

  !call calc_jet(nozzle, time)
  !if (dr_globalMe==MASTER_PE .and. mod(dr_nStep,20)==0) then
  !   write(*,'(a,2es11.3, f7.2)') '      (p, rho, M)=', &
  !   sim(nozzle)%pressure, sim(nozzle)%density, sim(nozzle)%mach
  !endif

  !write(*,*) 'blockCount:', blockCount
  nPtProc=0
  !allocate(pos(nPtProc,MDIM))
  !write(*,*) '[Heat] nPtProc:', nPtProc
  if (time.ge.sim(nozzle)%tOn .and. time.lt.(sim(nozzle)%tOn+sim(nozzle)%duration)) then

     call RANDOM_SEED(put=pt_randSeed)
     do blkInd=1,blockCount
        blockID = blockList(blkInd)
        call Grid_getDeltas(blockID,del)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        call Timers_start('Heat_fillnozzle')
        call Heat_fillnozzle(blockID,dt,time)
        call Timers_stop('Heat_fillnozzle')

     enddo
     ! Add shock and jet particles when the jet is on
     if (time.ge.(sim(nozzle)%tOn) .and. & 
         time.lt.(sim(nozzle)%tOn+0.99*sim(nozzle)%duration)) then
        gcMask = .false.
        gcMask(DENS_VAR) = .true.
        gcMask(PRES_VAR) = .true.
        gcMask(EINT_VAR) = .true.
        gcMask(ENER_VAR) = .true.
        gcMask(VELX_VAR:VELZ_VAR) = .true.
        gcMask(JET_SPEC) = .true.
        gcMask(ISM_SPEC) = .true.
        call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.true.,&
             maskSize=NUNK_VARS,mask=gcMask)

        ! Add shock particles
        ! nPtProc and pos were updated in Heat_fillnozzle
        call Timers_start('Particles_addNew')
        call Particles_addNew(nPtProc, pos, 1.0, addNewSuccess)
        call Timers_stop('Particles_addNew')
        !write(*,*) '[Heat2] nPtProc:', nPtProc
        ! Add new particles at the surfaces of the nozzle
        call Timers_start('Particles_addNew_nozzle')
        if (dr_globalMe==MASTER_PE) then

           call RANDOM_NUMBER(prob)
           !call MPI_Bcast(prob,1,MPI_DOUBLE_PRECISION,MASTER_PE,MPI_COMM_WORLD,ierr)

           !write(*,'(A6, 2f9.5)') 'prob', prob, 1.0/sim_ptAddPeriod*dt
           !write(*,'(i5, f9.5, es11.3)') pt_meshMe, prob, dtNew
           nPtNoz = int(1.0/sim_ptAddPeriod*dt)
           if (prob .le. 1.0/sim_ptAddPeriod*dt-nPtNoz) then
              nPtNoz = nPtNoz+1
           endif
           if (nPtNoz .gt. 0) then
              allocate(posNoz(nPtNoz,MDIM))
              call pt_getRandomPos(nPtNoz, posNoz)

              call Particles_addNew(nPtNoz, posNoz, 0.0, addNewSuccess)
              deallocate(posNoz)
           else
               !write(*,*) '[Heat] no pos', nPtProc, pos
               call Particles_addNew(0, pos, 0.0,  addNewSuccess)
               !write(*,*) '[Heat] no pos2'
           endif
        else
           !write(*,*) '[Heat] not masterpe', nPtProc, pos
           call Particles_addNew(0, pos, 0.0, addNewSuccess)
           !write(*,*) '[Heat] not masterpe2'
        endif
        call Timers_stop('Particles_addNew_nozzle')
        call RANDOM_SEED(get=pt_randSeed)
     endif

     !deallocate(pos)


     !write(*,'(i5, A28, i5)') pt_meshMe, 'After addNew, pt_numLocal=', pt_numLocal
     !write(*,'(i5, A28, i5)') pt_meshMe, 'After addNew, p_count    =', p_count
  endif

  return

end subroutine Heat


