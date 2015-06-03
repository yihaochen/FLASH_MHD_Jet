!!****if* source/Driver/DriverMain/Driver_computeDt
!!
!! NAME
!!
!!  Driver_computeDt
!!
!!
!! SYNOPSIS
!!
!!  Driver_computeDt(integer(IN) :: nbegin,
!!                  integer(IN) :: nstep, 
!!                  real(IN)    :: simTime,
!!                  real(IN)    :: dtOld,
!!                  real(OUT)   :: dtNew)
!!
!! DESCRIPTION
!!
!!  Determine the stability-limited time step.
!!  This timestep is determined using information from the included
!!  physics modules - many different timestep limiters are polled.
!!
!!  The global driver might use a different (hopefully smaller) time
!!  step, to match a file write time (tplot or trstr) or if the
!!  simulation end time has been reached; such possibilities are
!!  not considered here.
!!
!! ARGUMENTS
!!
!!  nbegin - first step of the simulation (nbegin is only used
!!              to determine if a label header should be written to
!!              the screen)
!!  nstep - current step of the simulation
!!  simTime - current simulation time of the run
!!  dtOld - the dt from the timestep that we just finished 
!!         (it's old because we be using dtOld to calculate 
!!          and return the dt for the next timestep (dtNew)
!!  dtNew - returned value of the dt calculated for the next timestep
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!! The calls to units currently not included in the code are commented out.
!!
!!
!!
!!*** 
#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_computeDt(nbegin, nstep, &
                    simTime, dtOld, dtNew)

  use Driver_data, ONLY : dr_dtMin,dr_dtMax, dr_tstepChangeFactor, &
                          dr_redshift, dr_useRedshift,             &
                          dr_printTStepLoc,                        &
                          dr_dtSTS, dr_useSTS, dr_globalMe, dr_globalComm,&
                          dr_dtAdvect, dr_dtDiffuse, dr_dtHeatExch
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getDeltas, &
    Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getSingleCellCoords
  use Hydro_interface, ONLY : Hydro_computeDt, Hydro_consolidateCFL
  use Stir_interface, ONLY: Stir_computeDt
  use Cosmology_interface, ONLY: Cosmology_computeDt
  use Cool_interface, ONLY : Cool_computeDt
  use Heat_interface, ONLY : Heat_computeDt
  use Heatexchange_interface, ONLY : Heatexchange_computeDt
  use Diffuse_interface, ONLY : Diffuse_computeDt 
  use Burn_interface, ONLY : Burn_computeDt
  use RadTrans_interface, ONLY: RadTrans_computeDt
  use Particles_interface, ONLY: Particles_computeDt

  use SolidMechanics_interface, ONLY : SolidMechanics_computeDt
  use IncompNS_interface, ONLY : IncompNS_computeDt
! <- ychen 06-2015
  use IO_data, ONLY : io_nextCheckPointStep
! ychen ->

  implicit none

#include "constants.h"
#include "Flash.h"
  include "Flash_mpi.h"

  integer, intent(IN) :: nbegin, nstep
  real,    intent(IN) :: simTime    !! current simulation time
  real, intent(IN) :: dtOld      !! last time step we used
  real, intent(OUT):: dtNew      !! the new timestep we get. to be returned.
 

  ! Local variables and functions
  integer :: i, error, blockID, numLeafBlocks, iout

  !! This is arbitrarily fixed. Users that add more units that influence the time
  !! should change this.

  integer, parameter :: nUnits = 15
  real, PARAMETER :: MAX_TSTEP = huge(1.0)

  
  real    :: dtModule(2,nUnits), dtLocal(2,nUnits)
  integer :: dtMinLoc(5), lminloc(5,nUnits), ngmin, pgmin
  integer :: status(MPI_Status_Size)

  logical :: gcell = .true.
  real, DIMENSION(MDIM) :: coords

  integer, dimension(MAXBLOCKS) :: blockList

  real, dimension(nUnits) :: tstepOutput
  character (len=20), DIMENSION(nUnits) :: &
                         limiterName, limiterNameOutput

  !!prepatory data structures for passing coords to timestep routines
  real, dimension(MDIM) :: del
  integer, dimension(MDIM) :: index

#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_ILO_GC:GRID_IHI_GC) :: xLeft,xRight,xCenter
  real,dimension(GRID_JLO_GC:GRID_JHI_GC) :: yLeft,yRight,yCenter
  real,dimension(GRID_KLO_GC:GRID_KHI_GC) :: zLeft,zRight,zCenter
  real, dimension(GRID_ILO_GC:GRID_IHI_GC) ::  dx, uxgrid
  real, dimension(GRID_JLO_GC:GRID_JHI_GC) ::  dy, uygrid
  real, dimension(GRID_KLO_GC:GRID_KHI_GC) ::  dz, uzgrid
#else
  real, allocatable,dimension(:)::&
       dx,uxgrid,dy,uygrid,dz,uzgrid
  real, allocatable,dimension(:)::xLeft,xRight,xCenter
  real, allocatable,dimension(:)::yLeft,yRight,yCenter
  real, allocatable,dimension(:)::zLeft,zRight,zCenter

#endif

  !arrays which hold the starting and ending indicies of a block
  integer,dimension(2,MDIM)::blkLimits,blkLimitsGC

  !!coordinate infomration to be passed into physics  
  real, pointer :: solnData(:,:,:,:)
  integer :: isize,jsize,ksize

  logical :: printTStepLoc
  integer :: itempLimit = 0
  integer, parameter :: HYDRO=1,BURN=2,GRAV=3,HEAT=4,COOL=5,TEMP=6,&
                        PART=7,DIFF=8,COSMO=9,STIR=10,HEATXCHG=11, &
                        RADTRANS=12,STS=13,INS=14,SOLID=15
  logical :: printToScrn
  real :: extraHydroInfo
  character (len=20) :: cflNumber
  real :: extraHydroInfoMin

  ! Initializing extraHydroInfo to zero:
  extraHydroInfo = 0.
  extraHydroInfoMin = 1.e10 !temporary large fake CFL for comparison


  data limiterName(HYDRO) /'dt_hydro'/
  data limiterName(HEAT) /'dt_Heat'/
  data limiterName(PART) /'dt_Part '/
  data limiterName(BURN) /'dt_Burn '/
  data limiterName(COOL) /'dt_Cool '/
  data limiterName(TEMP) /'dt_Temp '/
  data limiterName(DIFF) /'dt_Diff '/
  data limiterName(COSMO) /'dt_Cosm'/
  data limiterName(STIR) /'dt_Stir'/
  data limiterName(HEATXCHG) /'dt_HeatXchg'/
  data limiterName(RADTRANS) /'dt_RadTrans'/
  data limiterName(STS)  /'dt_STS'/
  data cflNumber  /'CFL'/


  !            Find the local minimum timestep among the included physics
  !            modules for locally stored blocks.
  
  !            Initialize all timestep variables.




  printTStepLoc = dr_printTStepLoc
  
  dtMinLoc(:) = 0
  lminloc(:,:) = 0
  lminloc(NDIM+1:MDIM,:) = 1

  do i = 1, nUnits
     dtLocal(1,i) = MAX_TSTEP
     dtLocal(2,i) = real(dr_globalMe)
  enddo

  !            Loop over all local leaf-node blocks

  call Hydro_consolidateCFL()
  
  call Grid_getListOfBlocks(LEAF,blockList, numLeafBlocks)
  
  do i = 1, numLeafBlocks
     !!Get the coordinate information for all the
     call Grid_getBlkIndexLimits(blockList(i),blkLimits,blkLimitsGC)
     isize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     jsize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     ksize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
#ifndef FIXEDBLOCKSIZE
     allocate(xLeft(isize))
     allocate(xRight(isize))
     allocate(xCenter(isize))
     allocate(dx(isize))
     allocate(uxgrid(isize))
     allocate(yLeft(jsize))
     allocate(yRight(jsize))
     allocate(yCenter(jsize))
     allocate(dy(jsize))
     allocate(uygrid(jsize))
     allocate(zLeft(ksize))
     allocate(zRight(ksize))
     allocate(zCenter(ksize))
     allocate(dz(ksize))
     allocate(uzgrid(ksize))
#endif
#ifdef DEBUG_DRIVER
     print*,'before calling get coordinates',isize,gcell
#endif
     call Grid_getCellCoords(IAXIS,blockList(i),CENTER,gcell,xCenter,isize)
     call Grid_getCellCoords(IAXIS,blockList(i),LEFT_EDGE,gcell,xLeft,isize)
     call Grid_getCellCoords(IAXIS,blockList(i),RIGHT_EDGE,gcell,xRight,isize)

#ifdef DEBUG_DRIVER
     print*,'before calling get coordinates',jsize,gcell
#endif
     if (NDIM > 1) then
        call Grid_getCellCoords(JAXIS,blockList(i),CENTER,gcell,yCenter,jsize)
        call Grid_getCellCoords(JAXIS,blockList(i),LEFT_EDGE,gcell,yLeft,jsize)
        call Grid_getCellCoords(JAXIS,blockList(i),RIGHT_EDGE,gcell,yRight,jsize)

        if (NDIM > 2) then
#ifdef DEBUG_DRIVER
           print*,'before calling get coordinates',ksize,gcell
#endif
           call Grid_getCellCoords(KAXIS,blockList(i),CENTER,gcell,zCenter,ksize)
           call Grid_getCellCoords(KAXIS,blockList(i),LEFT_EDGE,gcell,zLeft,ksize)
           call Grid_getCellCoords(KAXIS,blockList(i),RIGHT_EDGE,gcell,zRight,ksize)           
        endif
     endif

     call Grid_getDeltas(blockList(i), del)
     dx(:) = del(1)
     dy(:) = del(2)
     dz(:) = del(3)

     uxgrid(:) = 0
     uygrid(:) = 0
     uzgrid(:) = 0

     call Grid_getBlkPtr(blockList(i),solnData)

     ! hydro
#ifdef DEBUG_DRIVER
     print*,'going to call Hydro timestep'
#endif
     !extraHydroInfo = 0.
     call Hydro_computeDt (blockList(i), &
                           xCenter, dx, uxgrid, &
                           yCenter, dy, uygrid, &
                           zCenter, dz, uzgrid, &
                           blkLimits,blkLimitsGC,  &
                           solnData,      &
                           dtLocal(1,HYDRO), lminloc(:,HYDRO), &
                           extraInfo=extraHydroInfo )

     !! Extra CFL information
     if (extraHydroInfo .ne. 0.) then
        if (extraHydroInfo <= extraHydroInfoMin) then
           extraHydroInfoMin = extraHydroInfo
        endif
     else !if extraHydroInfo == 0.
        extraHydroInfoMin = extraHydroInfo
     endif

     call Stir_computeDt ( blockList(i),  &
                           blkLimits,blkLimitsGC,  &
                           solnData,               &
                          dtLocal(1,STIR), lminloc(:,STIR) )
   

#ifdef DEBUG_DRIVER
     print*,'returned from hydro timestep'
#endif

     call Burn_computeDt ( blockList(i),  &
                           blkLimits,blkLimitsGC,  &
                           solnData,               &
                           dtLocal(1,BURN), lminloc(:,BURN) )
     
!!$     call Gravity_computeDt ( blockList(i), dr_globalMe, &
!!$                           xCenter,xLeft,xRight, dx, uxgrid, &
!!$                           yCenter,yLeft,yRight, dy, uygrid, &
!!$                           zCenter,zLeft,zRight, dz, uzgrid, &
!!$                           blkLimits,blkLimitsGC,  &
!!$                           solnData,      &
!!$                           dtLocal(1,GRAV), lminloc(:,GRAV) )
     
     
     call Heat_computeDt ( blockList(i),  &
                           xCenter, dx, uxgrid, &
                           yCenter, dy, uygrid, &
                           zCenter, dz, uzgrid, &
                           blkLimits,blkLimitsGC,  &
                           solnData,      &
                           dtLocal(1,HEAT), lminloc(:,HEAT) )

     call Heatexchange_computeDt ( blockList(i),  &
                           blkLimits,blkLimitsGC,  &
                           solnData,      &
                           dtLocal(1,HEATXCHG), lminloc(:,HEATXCHG) )

     call RadTrans_computeDt(blockList(i), blkLimits,blkLimitsGC, &
          solnData, dtLocal(1,RADTRANS), lminloc(:,RADTRANS) )

     call Cool_computeDt ( blockList(i),  &
                           blkLimits,blkLimitsGC,  &
                           solnData,      &
                           dtLocal(1,COOL), lminloc(:,COOL) )

     call Particles_computeDt &
          ( blockList(i), dtLocal(1,PART), lminloc(:,PART))


     call Diffuse_computeDt ( blockList(i),  &
                           xCenter,xLeft,xRight, dx, uxgrid, &
                           yCenter,yLeft,yRight, dy, uygrid, &
                           zCenter,zLeft,zRight, dz, uzgrid, &
                           blkLimits,blkLimitsGC,  &
                           solnData,      &
                           dtLocal(1,DIFF), lminloc(:,DIFF) )
     
!!$     call Cosmo_timestep ( blockList(i),  &
!!$                           xCenter,xLeft,xRight, dx, uxgrid, &
!!$                           yCenter,yLeft,yRight, dy, uygrid, &
!!$                           zCenter,zLeft,zRight, dz, uzgrid, &
!!$                           blkLimits,blkLimitsGC,  &
!!$                           solnData,      &
!!$                           dtLocal(1,COSMO), lminloc(:,COSMO) )
      call Cosmology_computeDt(dtLocal(1,COSMO))


      

      !! Super time step
      if (dr_useSTS) then
         dtLocal(1,STS) = dr_dtSTS
      endif

#ifndef FIXEDBLOCKSIZE
     deallocate(xCenter)
     deallocate(xLeft)
     deallocate(xRight)
     deallocate(dx)
     deallocate(uxgrid)
     deallocate(yCenter)
     deallocate(yLeft)
     deallocate(yRight)
     deallocate(dy)
     deallocate(uygrid)
     deallocate(zCenter)
     deallocate(zLeft)
     deallocate(zRight)
     deallocate(dz)
     deallocate(uzgrid)
#endif
#ifdef DEBUG_DRIVER
     print*,'release blockpointer'
#endif

     call Grid_releaseBlkPtr(blockList(i),solnData)
  enddo

!!$  !! Choose the smallest CFL for screen output
  extraHydroInfo = 0.
  call MPI_AllReduce (extraHydroInfoMin, extraHydroInfo, 1, & 
       FLASH_REAL, MPI_MIN, dr_globalComm, error)


  ! IncompNS:
  call IncompNS_computeDt(dtLocal(1,INS),lminloc(:,INS))


  ! SolidMechanics:
  call SolidMechanics_computeDt(dtLocal(1,SOLID))


  ! DEV: we disabled temperature timestep limiter for now.  
  ! The old temperature was not updated with the refinement, 
  ! so dT/T was precomputed and the number of blocks may not be
  ! the same as there are now.
!!$  if (itempLimit == 1) then
!!$     do blockID = 1,  MAXBLOCKS
!!$       call Driver_computeDtTemp(dr_globalMe, dtOld, dtLocal(1,6), &
!!$             lminloc(1,6), blockID)
!!$     enddo
!!$  endif
  


  !            Find the minimum timestep across all processors and all
  !            modules.

  call MPI_AllReduce (dtLocal(1,1), dtModule(1,1), nUnits, & 
       MPI_2Double_Precision, MPI_MinLoc, dr_globalComm, error)

  dtNew    = huge(dtNew)    ! dt will hold the minimum timestep
  ngmin = 1                 ! ngmin will hold the winning module #
  pgmin = MASTER_PE         ! pgmin will hold the winning PE #
  
!!$  do i = 1, nUnits-1
!!$     if (dtModule(1,i) < dtNew) then
!!$        dtNew = dtModule(1,i)
!!$        pgmin = dtModule(2,i)
!!$        ngmin = i
!!$     endif
!!$  enddo

  do i = 1, nUnits
     if ((i .ne. STS) .and. (i .ne. DIFF)) then
        if (dtModule(1,i) < dtNew) then
           dtNew = dtModule(1,i)
           pgmin = dtModule(2,i)
           ngmin = i
        endif
     endif
  enddo

  ! Save it to hydro's advection time scale
  ! Note: This dr_dtAdvect is the minimum timestep from all physics units,
  !       e.g., Hydro, Stir, Burn, Heat, Cool, Particle, and Cosmology,
  !       except from DIFF and STS.
  dr_dtAdvect = dtNew

  ! Do it one more time
  if (dtModule(1,DIFF) < dtNew) then
     dtNew = dtModule(1,DIFF)
     pgmin = dtModule(2,DIFF)
     ngmin = DIFF
  endif

  ! Save it to hydro's diffusion time scale
  dr_dtDiffuse = dtModule(1,DIFF)
      
  ! have the processor that is determining the timestep limit broadcast the
  ! proc number, block number, and i,j,k of the zone that set the timestep
  ! to all processors

  dtMinLoc(:) = lminloc(:,ngmin)

  call MPI_Bcast(dtMinLoc(1), 5, MPI_INTEGER, pgmin, dr_globalComm, error)

  ! limit the timestep to increase by at most a factor of dr_tstepChangeFactor

  dtNew = min( dtNew, dtOld*dr_tstepChangeFactor )
  
  if (nstep .GE. 50) then   !! This is where Cellular starts to fail
!     print *, 'nstep = ',nstep
  endif


! <- ychen 06-2015
  ! Write a checkpoint file if dtNew is smaller than dr_dtMin (dtMin set in flash.par)
  if (dtNew .LT. dr_dtMin) then
     io_nextCheckpointStep = nStep+1
  endif
! ychen ->

  ! Use dr_dtmin and dr_dtmax to limit the timestep.  If this makes the code
  ! unstable, it's not our fault.
  dtNew = min( max( dtNew, dr_dtMin ), dr_dtMax )


  if (printTStepLoc) then
         
     ! convert the dtMinLoc array into a physical position (x,y,z) where the
     ! timestep is being set.  dtMinLoc(5) is the processor number, dtMinLoc(4)
     ! is the blockID on that proc.
     coords(:) = 0.0

     if (dr_globalMe == dtMinLoc(5)) then

        if (dtMinLoc(4) > 0) then
           index(:)=dtMinLoc(1:MDIM)
           call Grid_getSingleCellCoords(index,dtMinLoc(4),CENTER, EXTERIOR,coords)
        else
           coords(:) = 999.0
        end if

        ! send this to the master processor
        if (dr_globalMe /= MASTER_PE) then

           call MPI_Send (coords(1), 3, FLASH_REAL, MASTER_PE, & 
                0, dr_globalComm, error)
           
        endif
        
     elseif (dr_globalMe == MASTER_PE) then
        
        call MPI_Recv (coords(1), 3, FLASH_REAL, dtMinLoc(5), 0, & 
             dr_globalComm, status, error)            
        
     endif
     
  endif
  
  ! Print out the time and next timestep.
  
  ! only print out the timestep from the limiters that are active
  iout = 0
  do i = 1, nUnits
     if (dtModule(1,i) /= MAX_TSTEP) then
        iout = iout + 1
        tstepOutput(iout) = dtModule(1,i)
        limiterNameOutput(iout) = limiterName(i)
     endif
  enddo

!!$print*,iout,nUnits;pause


  printToScrn = .true.
  if (printToScrn) then
  if (dr_globalMe == MASTER_PE) then


  if (extraHydroInfo .eq. 0.) then

     if (printTStepLoc) then
        
        if (nstep == nbegin) then

              if (.not. dr_useRedshift) then
                 write (*,803) 'n', 't', 'dt', 'x', 'y', 'z', &
                      (limiterNameOutput(i),i=1,iout)
              else
                 write (*,804) 'n', 't', 'z', 'dt', 'x', 'y', 'z', &
                      (limiterNameOutput(i),i=1,iout)
              endif
              
           endif
        
           if (.not. dr_useRedshift) then
              if (.not. dr_useSTS) then
                 write(*,801) nstep, simTime, dtNew, coords(1), coords(2), &
                      coords(3), (tstepOutput(i),i=1,iout)
              else
                 write(*,801) nstep, simTime, max(dtNew,dr_dtSTS), coords(1), coords(2), &
                      coords(3), (tstepOutput(i),i=1,iout)
              endif

           else
              if (.not. dr_useSTS) then
                 write(*,802) nstep, simTime, dr_redshift, dtNew, coords(1), &
                      coords(2), coords(3), (tstepOutput(i),i=1,iout)
              else
                 write(*,802) nstep, simTime, dr_redshift, max(dtNew,dr_dtSTS), coords(1), &
                      coords(2), coords(3), (tstepOutput(i),i=1,iout)
              endif
           endif
           
        else
        
           if (nstep .eq. nbegin) then
           
              if (.not. dr_useRedshift) then
                 write (*,903) 'n', 't', 'dt', (limiterNameOutput(i),i=1,iout)
              else
                 write (*,904) 'n', 't', 'z', 'dt', (limiterNameOutput(i),i=1,iout)
              endif
           
           endif
        
           if (.not. dr_useRedshift) then
              if (.not. dr_useSTS) then
                 write(*,901) nstep, simTime, dtNew, (tstepOutput(i),i=1,iout)
              else
                 write(*,901) nstep, simTime, max(dtNew,dr_dtSTS), (tstepOutput(i),i=1,iout)
              endif
           else
              if (.not. dr_useSTS) then
                 write(*,902) nstep, simTime, dr_redshift, dtNew, &
                      (tstepOutput(i),i=1,iout)
              else
                 write(*,902) nstep, simTime, dr_redshift, max(dtNew,dr_dtSTS), &
                      (tstepOutput(i),i=1,iout)
              endif

           endif

        endif

     else ! elseif (extraHydroInfo .ne. 0.) then

        if (printTStepLoc) then
           if (nstep == nbegin) then

              if (.not. dr_useRedshift) then
                 write (*,803) 'n', 't', 'dt', 'x', 'y', 'z', &
                      (limiterNameOutput(i),i=1,iout),cflNumber
              else
                 write (*,804) 'n', 't', 'z', 'dt', 'x', 'y', 'z', &
                      (limiterNameOutput(i),i=1,iout),cflNumber
              endif
              
           endif
        
           if (.not. dr_useRedshift) then
              if (.not. dr_useSTS) then
                 write(*,801) nstep, simTime, dtNew, coords(1), coords(2), &
                      coords(3), (tstepOutput(i),i=1,iout), extraHydroInfo
              else
                 write(*,801) nstep, simTime, max(dtNew,dr_dtSTS), coords(1), coords(2), &
                      coords(3), (tstepOutput(i),i=1,iout), extraHydroInfo
              endif

           else
              if (.not. dr_useSTS) then
                 write(*,802) nstep, simTime, dr_redshift, dtNew, coords(1), &
                      coords(2), coords(3), (tstepOutput(i),i=1,iout), extraHydroInfo
              else
                 write(*,802) nstep, simTime, dr_redshift, max(dtNew,dr_dtSTS), coords(1), &
                      coords(2), coords(3), (tstepOutput(i),i=1,iout), extraHydroInfo
              endif
           endif
           
        else
        
           if (nstep .eq. nbegin) then
           
              if (.not. dr_useRedshift) then
                 write (*,903) 'n', 't', 'dt', (limiterNameOutput(i),i=1,iout),cflNumber
              else
                 write (*,904) 'n', 't', 'z', 'dt', (limiterNameOutput(i),i=1,iout),cflNumber
              endif
           
           endif
        
           if (.not. dr_useRedshift) then
              if (.not. dr_useSTS) then
                 write(*,901) nstep, simTime, dtNew, (tstepOutput(i),i=1,iout), extraHydroInfo
              else
                 write(*,901) nstep, simTime, max(dtNew,dr_dtSTS), (tstepOutput(i),i=1,iout), extraHydroInfo
              endif
           else
              if (.not. dr_useSTS) then
                 write(*,902) nstep, simTime, dr_redshift, dtNew, &
                      (tstepOutput(i),i=1,iout), extraHydroInfo
              else
                 write(*,902) nstep, simTime, dr_redshift, max(dtNew,dr_dtSTS), &
                      (tstepOutput(i),i=1,iout), extraHydroInfo
              endif

           endif

        endif

     endif


  endif

endif ! end of printToScrn
      

801 format (1X, I7, 1x, ES10.4, 1x, ES10.4, 2x, '(', ES10.3, ', ', &
            ES 10.3, ', ', ES10.3, ')', ' | ', 11(1X, :, ES9.3),1x,ES10.3)
802 format (1X, I7, 1x, ES10.4, 1x, F8.3, 1x, ES10.4, 2x, '(', ES9.3, ', ', &
            ES 9.3, ', ', ES9.3, ')', ' | ', 11(1X, :, ES9.3),1x,ES10.3)
803 format (1X, A7, 1x, A10, 1x, A10, 2x, '(', A10, ', ', A10, ', ', A10, ')', &
            ' | ', 11(1X, :, A9),1x,A10)
804 format (1X, A7, 1x, A10, 1x, A7, 1x, A10, 2x, '(', A9, ', ', A9, ', ', &
         A9, ')', ' | ', 11(1X, :, A9),1x,A10)

901 format (1X, I7, 1X, ES10.4, 1x, ES10.4, ' | ', 11(1X, :, ES11.5),1x,ES10.4)
902 format (1X, I7, 1X, ES10.4, 1x, F8.3, 1x, ES10.4, ' | ', 11(1X, :, ES11.5),1x,ES10.4)
903 format (1X, A7, 1x, A10   , 1x, A10,    ' | ', 11(1X, :, A11),1x,A10)
904 format (1X, A7, 1x, A10   , 1x, A7, 1x, A10,    ' | ', 11(1X, :, A11),1x,A10)

  return
end subroutine Driver_computeDt
