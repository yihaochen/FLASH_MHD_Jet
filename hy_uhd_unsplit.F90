!!****if* source/physics/Hydro/HydroMain/unsplit/MHD_StaggeredMesh/hy_uhd_unsplit
!!
!! NAME
!!
!!  hy_uhd_unsplit
!!
!! SYNOPSIS
!!
!!  call hy_uhd_unsplit(integer (IN) :: blockCount,
!!                  integer (IN) :: blockList(blockCount),
!!                  real    (IN) :: dt,
!!                  real    (IN) :: dtOld  )
!!
!! DESCRIPTION
!!
!!  Performs MHD update in a directionally unsplit fashion over a set
!!  of blocks. Divergence cleaning control for the magnetic fields is
!!  handled by using the staggered mesh algorithm with a call to
!!  hy_uhd_staggeredDivb. Before calling this routine, electric fields are
!!  to be calculated with a call to hy_uhd_getElectricFields.
!!  blockList is an integer array of size blockCount that contains 
!!  the blocks over which to update.
!!  dt gives the timestep over which to advance, and timeEndAdv gives the
!!  simulation time at the end of the update.
!! 
!!  This routine performs a guardcell fill and for each block: 
!!   - Eos_wrapped call is applied to the guard cells where necessary;
!!   - computes fluxes using a call to hy_uhd_getFaceFlux
!!   - if we're not doing flux correction (as controlled by the flux_correct
!!     runtime parameter), then we update all the cell values from the fluxes 
!!     (with a call to hy_uhd_unsplitUpdate), otherwise, we update just cells 
!!     not on the boundaries, and save fluxes for cells on the boundary;
!!   - and finally, we apply an eos to the block.
!! 
!!  After the main block loop, if doing flux correction, we have
!!  the Grid correct boundary fluxes for all blocks where approriate,
!!  and do another loop over blocks, updating the cell values for
!!  cells on the block boundaries using the corrected fluxes, and
!!  apply an eos on the block. 
!!  The same is true for the electric fields correction where
!!  there are block boundaries that are sharing different levels of
!!  refinements.
!!
!!
!! REFERENCES
!!
!!  * Lee, D., Ph.D. Dissertation, Univ. of MD, 2006
!!  * Lee, D. and Deane, A., "An Unsplit Staggered Mesh Scheme for Multidimensional
!!                            Magnetohydrodynamics", 228 (2009), 952-975, JCP
!!  * Lee, D., "A solution accurate, efficient and stable unsplit staggered mesh scheme 
!!              for three dimensional magnetohydrodynamics", 243 (2013), 269-292, JCP
!!
!!
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList  -  array holding local IDs of blocks on which to advance
!!  dt         -  timestep
!!  dtOld      -  old timestep (needed by some implementations for temporal extrapolations
!!                              of gravity)
!!
!!***

!!REORDER(4): U, scrch_Ptr, scrchFace[XYZ]Ptr, fl[xyz]

#ifdef DEBUG_ALL
#define DEBUG_MHD
#endif
#define DEBUG_GRID_GCMASK

Subroutine hy_uhd_unsplit ( blockCount, blockList, dt, dtOld )

  use Hydro_data, ONLY : hy_fluxCorrect,      &
                         hy_gref,             &
                         hy_useGravity,       &
                         hy_order,            &
                         hy_units,            &
                         hy_gcMaskSize,       &
                         hy_gcMask,           &
                         hy_unsplitEosMode,   &
                         hy_eosModeGc,        &
                         hy_eosModeAfter,     &
                         hy_updateHydroFluxes,&
                         hy_geometry,         &
                         hy_fluxCorVars,      &
                         hy_cfl,              &
                         hy_cfl_original,     &
                         hy_numXN,            &
                         hy_fullRiemannStateArrays,    &
                         hy_fullSpecMsFluxHandling,   &
                         hy_dtmin,            &
                         hy_simTime,          &
                         hy_simGeneration,    &
                         hy_shockDetectOn,    &
                         hy_doUnsplitLoop0,   &
                         hy_killdivb,         &
                         hy_killdivb8w,       &
                         hy_forceHydroLimit,  &
                         hy_useBiermann,      &
                         hy_useBiermann1T,    &
                         hy_biermannSource

  use Driver_interface, ONLY : Driver_abortFlash, Driver_getSimTime
! <- ychen 10-2014
  use Driver_data, ONLY : dr_simTime, dr_globalMe, dr_nStep
  use Simulation_jetNozzleUpdate, ONLY : sim_jetNozzleUpdate
  use Simulation_data
! ychen ->

  use hy_uhd_interface, ONLY : hy_uhd_getRiemannState,  &
                               hy_uhd_getFaceFlux,      &
                               hy_uhd_unsplitUpdate,    &
                               hy_uhd_unitConvert,      &
                               hy_uhd_energyFix,        &
                               hy_uhd_prepareNewGravityAccel,&
                               hy_uhd_putGravityUnsplit,&
                               hy_uhd_addGravityUnsplit,&
                               hy_uhd_shockDetect,      &
                               hy_uhd_biermannSource,   &
                               hy_uhd_getElectricFields,&
                               hy_uhd_staggeredDivb,    &
                               hy_uhd_unsplitUpdateMultiTemp,&
                               hy_uhd_multiTempAfter
  use hy_memInterface, ONLY :  hy_memAllocScratch,      &
                               hy_memDeallocScratch,    &
                               hy_memGetBlkPtr,         &
                               hy_memReleaseBlkPtr


  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes,    &
                             Grid_conserveField,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkData,        &
                             Grid_getBlkNeighLevels

  use Eos_interface, ONLY : Eos_wrapped

  use Logfile_interface, ONLY : Logfile_stampVarMask, Logfile_stamp

  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"
#include "UHD.h"


  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList
  real,    INTENT(IN) :: dt, dtOld
  !! -----------------------------------------------------

  integer, dimension(MDIM) :: datasize
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, parameter :: level=0
  integer :: i,blockID
  integer :: ix,iy,iz
  real, dimension(MDIM) :: del
  logical :: gcMask(hy_gcMaskSize)
  logical :: halfTimeAdvance

#ifdef FIXEDBLOCKSIZE
  real :: flx(NFLUXES,&
              GRID_ILO_GC:GRID_IHI_GC,     &
              GRID_JLO_GC:GRID_JHI_GC,     &
              GRID_KLO_GC:GRID_KHI_GC)
  real :: fly(NFLUXES,&
              GRID_ILO_GC:GRID_IHI_GC,     &
              GRID_JLO_GC:GRID_JHI_GC,     &
              GRID_KLO_GC:GRID_KHI_GC)
  real :: flz(NFLUXES,&
              GRID_ILO_GC:GRID_IHI_GC,     &
              GRID_JLO_GC:GRID_JHI_GC,     &
              GRID_KLO_GC:GRID_KHI_GC)

  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: &
       gravX,gravY,gravZ
  real :: faceAreas(GRID_ILO_GC:GRID_IHI_GC,     &
                    GRID_JLO_GC:GRID_JHI_GC,     &
                    GRID_KLO_GC:GRID_KHI_GC)
#else
  real, allocatable, dimension(:,:,:,:) :: flx,fly,flz
  real, allocatable, dimension(:,:,:) :: gravX, gravY, gravZ
  real, allocatable :: faceAreas(:,:,:)
#endif

  real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr
  real, pointer, dimension(:,:,:,:) :: scrch_Ptr
  real, pointer, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
  integer :: firstScrCtr, lastScrCtr
#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif

  real, pointer, dimension(:,:,:,:) :: U

  integer :: updateMode ! will be set to one of UPDATE_ALL, UPDATE_INTERIOR, UPDATE_BOUND, etc.
  logical,dimension(MAXBLOCKS) :: blockNeedsFluxCorrect, blockMustStoreFluxes
  logical, parameter :: updateEarly=.FALSE. !update blocks in the first loop for which it is possible?
  integer :: neighLev(-1:1, -K2D:K2D , -K3D:K3D)
  integer :: myRefine
! <- ychen 10-2014
  integer :: nozzle=1
! ychen ->
  !! End of data declaration ***********************************************

  call Timers_start("Head")

#ifdef FLASH_GRID_PARAMESH2
  call Driver_abortFlash("The staggeredMesh MHD solver only works with PARAMESH 3 or 4!")
#endif



#ifdef FLASH_GRID_UG
  hy_fluxCorrect = .false.
#endif
  if (hy_killdivb) then
     if (hy_killdivb8w) then
        call Driver_abortFlash("You have both CT and 8-wave on. Please select one as true.")
     endif 
  endif 
  !! ***************************************************************************
  !! Shock detection before MHD                                                *
  !! ***************************************************************************
  !! Call shock detect algorithm to determine tagging shocks before hydro begins;
  !! other preparations of UNK data if necessary.
  if (hy_shockDetectOn) then

     !! Call guardcell filling to properly detect shocks
     gcMask = .false.
     gcMask(DENS_VAR) = .true.
     gcMask(PRES_VAR) = .true.
     gcMask(GAMC_VAR) = .true.
     gcMask(VELX_VAR:VELZ_VAR) = .true.
#ifdef CFL_VAR
     gcMask(CFL_VAR)  = .true.
#endif
#if NSPECIES > 1
     gcMask(SPECIES_BEGIN:SPECIES_END) = .true.
#endif

#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stampVarMask(gcMask, .FALSE., '[hy_uhd_unsplit]', 'gcWant[Detect]')
     end if
#endif

     call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.false.,&
          maskSize=NUNK_VARS,mask=gcMask,makeMaskConsistent=.false.,&
          doLogMask=.NOT.gcMaskLogged)

     hy_cfl = hy_cfl_original
  end if

  if (hy_doUnsplitLoop0) then

     do i=1,blockCount          !LOOP 0
        blockID = blockList(i)

        !! Detect shocks
        if (hy_shockDetectOn) call hy_uhd_shockDetect(blockID)

        if ( hy_units .NE. "NONE" .and. hy_units .NE. "none" ) then
           call hy_uhd_unitConvert(blockID,FWDCONVERT)
        endif

#if defined(GPRO_VAR)||defined(VOLX_VAR)||defined(VOLY_VAR)||defined(VOLZ_VAR)||defined(CFL_VAR)
        if (hy_updateHydroFluxes) then
           call Grid_getBlkPtr(blockID,U,CENTER)
#ifdef GPRO_VAR
           ! A tagging variable for Gaussian Process (GP) method.
           U(GPRO_VAR,:,:,:) = 0.
#endif
           !! -----------------------------------------------------------------------!
           !! Save old velocities ---------------------------------------------------!
           !! -----------------------------------------------------------------------!
#ifdef VOLX_VAR
           U(VOLX_VAR,:,:,:) = U(VELX_VAR,:,:,:)
#endif
#ifdef VOLY_VAR
           U(VOLY_VAR,:,:,:) = U(VELY_VAR,:,:,:)
#endif
#ifdef VOLZ_VAR
           U(VOLZ_VAR,:,:,:) = U(VELZ_VAR,:,:,:)
#endif
#ifdef CFL_VAR
           where (1.2*U(CFL_VAR,:,:,:) < hy_cfl_original .AND. U(CFL_VAR,:,:,:) .GE. 2.5)
              !! Slow recover (of factor of 1.2) to the original CFL once it gets to
              !! reduced to a smaller one in the presence of strong shocks.
              !! This variable CFL takes place in the following three cases using:
              !! (1) use_hybridOrder = .true.,
              !! (2) use_hybridOrder = .true., or
              !! (3) BDRY_VAR is defined and used for stationary objects.
              U(CFL_VAR,:,:,:) = 1.2*U(CFL_VAR,:,:,:)
           elsewhere (U(CFL_VAR,:,:,:) .GE. 2.5)
              U(CFL_VAR,:,:,:) = hy_cfl_original
           end where
#endif
           call Grid_releaseBlkPtr(blockID,U,CENTER)
        end if
#endif
     enddo
  endif


  !! ***************************************************************************
  !! Call guardcell filling with Eos before MHD                                *
  !! ***************************************************************************
#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     call Logfile_stampVarMask(hy_gcMask, .TRUE., '[hy_uhd_unsplit]', 'gcNeed')
  end if
#endif

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,doEos=.true.,eosMode=hy_eosModeGc,&
       maskSize=hy_gcMaskSize,mask=hy_gcMask,makeMaskConsistent=.true.,&
       doLogMask=.NOT.gcMaskLogged)


  !! ***************************************************************************
  !! First part of advancement                                                 *
  !! ***************************************************************************
  !! Loop over the blocks

  !! Retain the original cfl that may have been changed in some leaf blocks.
  if (hy_updateHydroFluxes) then
     if (1.2*hy_cfl < hy_cfl_original) then
        !! Slow recover (of factor of 1.2) to the original CFL once it gets to
        !! reduced to a smaller one in the presence of strong shocks.
        !! This variable CFL takes place in the following three cases using:
        !! (1) use_hybridOrder = .true.,
        !! (2) use_hybridOrder = .true., or
        !! (3) BDRY_VAR is defined and used for stationary objects.
        hy_cfl = 1.2*hy_cfl
     else
        hy_cfl = hy_cfl_original
     endif
     hy_dtmin = huge(1.0)
  endif

! <- ychen 10-2014
  call Timers_start('sim_jetNozzleUpdate')
  call sim_jetNozzleUpdate(nozzle, dr_simTime, dt)
  call Timers_stop('sim_jetNozzleUpdate')
  if (dr_globalMe==MASTER_PE .and. mod(dr_nStep,20)==0) then
     write(*,'(a,2es11.3, f7.2)') '      (p, rho, M)=', &
     sim(nozzle)%pressure, sim(nozzle)%density,&
     sim(nozzle)%velocity/sqrt(sim(nozzle)%gamma*sim(nozzle)%pressure/sim(nozzle)%density)
     write(*,'(a, 2es11.3)') '      (bz, bphi)=', sim(nozzle)%bz, sim(nozzle)%bphi
  endif
! ychen ->

#if defined(FLASH_UHD_3T)
  firstScrCtr = min(HY_VAR1_SCRATCHCTR_VAR,HY_XN01_SCRATCHCTR_VAR)
  lastScrCtr  = max(HY_VAR2_SCRATCHCTR_VAR,HY_XN06_SCRATCHCTR_VAR)
#else
  firstScrCtr = HY_VAR1_SCRATCHCTR_VAR
  lastScrCtr  = HY_VAR2_SCRATCHCTR_VAR
#endif
  call hy_memAllocScratch(SCRATCH_CTR,firstScrCtr,lastScrCtr-firstScrCtr+1, 0,0,0, &
       blockList(1:blockCount) )

  if (hy_fullRiemannStateArrays) then
     call hy_memAllocScratch(SCRATCH_FACEX,&
       min(HY_P01_FACEXPTR_VAR,HY_N01_FACEXPTR_VAR),&
       2*HY_SCRATCH_NUM, &
       0,1,2, &
       blockList(1:blockCount) )
#if NDIM > 1
     call hy_memAllocScratch(SCRATCH_FACEY,&
       min(HY_P01_FACEYPTR_VAR,HY_N01_FACEYPTR_VAR),&
       2*HY_SCRATCH_NUM,&
       0,1,2, &
       blockList(1:blockCount) )
#endif
#if NDIM > 2
     call hy_memAllocScratch(SCRATCH_FACEZ,&
       min(HY_P01_FACEZPTR_VAR,HY_N01_FACEZPTR_VAR),&
       2*HY_SCRATCH_NUM,&
       0,1,2, &
       blockList(1:blockCount) )
#endif
  end if

#if (NSPECIES+NMASS_SCALARS) > 0
     if (hy_fullSpecMsFluxHandling .AND. hy_fullRiemannStateArrays &
          .AND. NFACE_VARS > 0 .AND. NDIM > 1) then
        ! for hy_spcR:
        call hy_memAllocScratch(CENTER,&
             1,&
             hy_numXN,&
             2,0,0, &
             blockList(1:blockCount), &
             highSize=NDIM)
        ! for hy_spcL:
        call hy_memAllocScratch(SCRATCH_CTR,&
             1,&
             hy_numXN,&
             2,0,0, &
             blockList(1:blockCount), &
             highSize=NDIM)
     end if
#endif

  call Timers_stop("Head")

  do i=1,blockCount             !LOOP 1

     blockID = blockList(i)

     if (hy_fluxCorrect .AND. updateEarly) then
        ! Test whether neighbors are at different refinement levels, and if so,
        ! decide on what we have to do with this block for flux correction.
        call Grid_getBlkNeighLevels(blockID,neighLev)
        myRefine = neighLev(0,0,0)
        blockNeedsFluxCorrect(blockID) = (neighLev(-1,0,0) < myRefine .OR. neighLev(1,0,0) < myRefine)
        blockMustStoreFluxes (blockID) = (neighLev(-1,0,0) .NE. myRefine .OR. neighLev(1,0,0) .NE. myRefine)
#if NDIM > 1
        blockNeedsFluxCorrect(blockID) = blockNeedsFluxCorrect(blockID) &
             .OR. (neighLev(0,-1,0) < myRefine .OR. neighLev(0,1,0) < myRefine)
        blockMustStoreFluxes (blockID) = blockMustStoreFluxes(blockID) &
             .OR. (neighLev(0,-1,0) .NE. myRefine .OR. neighLev(0,1,0) .NE. myRefine)
#endif
#if NDIM > 2
        blockNeedsFluxCorrect(blockID) = blockNeedsFluxCorrect(blockID) &
             .OR. (neighLev(0,0,-1) < myRefine .OR. neighLev(0,0,1) < myRefine)
        blockMustStoreFluxes (blockID) = blockMustStoreFluxes(blockID) &
             .OR. (neighLev(0,0,-1) .NE. myRefine .OR. neighLev(0,0,1) .NE. myRefine)
#endif
     else if (hy_fluxCorrect) then
        blockNeedsFluxCorrect(blockID) = .TRUE.
        blockMustStoreFluxes (blockID) = .TRUE.
     else
        blockNeedsFluxCorrect(blockID) = .FALSE.
        blockMustStoreFluxes (blockID) = .FALSE.
     end if

     call Grid_getDeltas(blockID,del)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

     if (hy_fullRiemannStateArrays) then
        call hy_memGetBlkPtr(blockID,scrchFaceXPtr,SCRATCH_FACEX)
        if (NDIM > 1) call hy_memGetBlkPtr(blockID,scrchFaceYPtr,SCRATCH_FACEY)
        if (NDIM > 2) call hy_memGetBlkPtr(blockID,scrchFaceZPtr,SCRATCH_FACEZ)
     else
        allocate(scrchFaceXPtr(HY_NSCRATCH_VARS,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(scrchFaceYPtr(HY_NSCRATCH_VARS,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(scrchFaceZPtr(HY_NSCRATCH_VARS,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     endif

#if (NSPECIES+NMASS_SCALARS) > 0
     if (hy_fullSpecMsFluxHandling) then
        if (hy_fullRiemannStateArrays .AND. NFACE_VARS > 0 .AND.  NDIM > 1) then
           call hy_memGetBlkPtr(blockID,hy_spcR,CENTER)
           call hy_memGetBlkPtr(blockID,hy_spcL,SCRATCH_CTR)
        else
           allocate(  hy_SpcR(HY_NSPEC,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS),NDIM))
           allocate(  hy_SpcL(HY_NSPEC,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS),NDIM))
        end if
        allocate(hy_SpcSig(HY_NSPEC,blkLimits(LOW,IAXIS)-3:blkLimits(HIGH,IAXIS)+3,dataSize(JAXIS),dataSize(KAXIS),NDIM))
     end if
#endif

#ifndef FIXEDBLOCKSIZE
     allocate(flx(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(fly(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(flz(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(    gravX(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(    gravY(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(    gravZ(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(  faceAreas(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif

     !! ************************************************************************
     !! Get gravity
     gravX = 0.
     gravY = 0.
     gravZ = 0.
     if (hy_useGravity) then
        call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
        gravX = gravX/hy_gref
        gravY = gravY/hy_gref
        gravZ = gravZ/hy_gref
     endif


     if (hy_updateHydroFluxes) then
        !! ************************************************************************
        !! Calculate Riemann (interface) states
        !! Note: gravX(:,:,:) - gravity at n

#if (NSPECIES+NMASS_SCALARS) > 0
        if (hy_fullSpecMsFluxHandling) then
           hy_SpcL=0.
           hy_SpcR=0.
           hy_SpcSig=0.
        end if
#endif
        call Timers_start("getRiemannState")
        call hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,dt,del, &
                                    gravX(:,:,:),gravY(:,:,:),gravZ(:,:,:),&
                                    scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,&
                                    hy_SpcR,hy_SpcL,hy_SpcSig,&
                                    normalFieldUpdate=.false.)
        call Timers_stop("getRiemannState")

        !! DEV: DL-This note seems to be outdated and wrong for the optimized code.
        ! Note: Two different ways of handling gravity:
        ! 1. With gravity calculated potential at n+1/2, Riemann states do not include gravity
        !    source terms at this point, and will include them in hy_uhd_addGravityUnsplit later
        !    to the primitive Riemann variables (not available for conservative formulation).
        ! 2. With gravity extrapolated from n-1 & n states, gravity source terms have been
        !    included to Riemann states in conservative formulation in hy_uhd_getRiemannState.

     endif !! End of if (hy_updateHydroFluxes) then
#if (NSPECIES+NMASS_SCALARS) > 0
     if (hy_fullSpecMsFluxHandling) then
        deallocate(hy_SpcSig)
     end if
#endif


     !!************************************************************************************!!
     !!*** Begin multidimensional MHD reconstruction (Riemann state) calculations *********!!
     !!************************************************************************************!!
     !! Note:
     !! Don't do this part when pure hydro (NFACE_VARS = 0) is considered and when NDIM = 1
     !! Do this intermediate part for multidimensional problem only.
     !! One dimensional problem doesn't need to do this part as 1D setup
     !! doesn't need to use magnetic fields (facevars) at cell face centers.
     !!
#if NFACE_VARS > 0
#if NDIM > 1

     !! ************************************************************************
     !! Evolve facevars n+1/2 time step to get 2nd order accuracy
     flx = 0.
     fly = 0.
     flz = 0.
     if ((.not. hy_forceHydroLimit) .and. hy_order > 1) then
        call hy_memGetBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)
        !! Intermediate Godunov fluxes
        call Timers_start("getFaceFlux")
        call hy_uhd_getFaceFlux(blockID,blkLimits,blkLimitsGC,datasize,del,flx,fly,flz,&
                                scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,scrch_Ptr,&
                                hy_SpcR,hy_SpcL,&
                                lastCall=.FALSE.)
        call Timers_stop("getFaceFlux")
        call hy_memReleaseBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)
        if (hy_killdivb) then
           !! DivB=0 calls
           halfTimeAdvance = .true.
           !! Calculate Godunov fluxes one more time to calculate electric fields
           call hy_uhd_getElectricFields(blockID,blkLimits,blkLimitsGC,del,flx,fly,flz)
        endif
     endif
     !! End of evolution of face-centered magnetic fields by n+1/2 time step
     !! ************************************************************************

#ifndef FIXEDBLOCKSIZE
     deallocate(flx)
     deallocate(fly)
     deallocate(flz)
     deallocate(gravX)
     deallocate(gravY)
     deallocate(gravZ)
     deallocate(faceAreas)
#endif

     call hy_memReleaseBlkPtr(blockID,scrchFaceXPtr,SCRATCH_FACEX)
     if (NDIM > 1) call hy_memReleaseBlkPtr(blockID,scrchFaceYPtr,SCRATCH_FACEY)
     if (NDIM > 2) call hy_memReleaseBlkPtr(blockID,scrchFaceZPtr,SCRATCH_FACEZ)

#if (NSPECIES+NMASS_SCALARS) > 0
     if (hy_fullSpecMsFluxHandling) then
        call hy_memReleaseBlkPtr(blockID,hy_spcR,CENTER)
        call hy_memReleaseBlkPtr(blockID,hy_spcL,SCRATCH_CTR)
     end if
#endif

  enddo
  !! End of the first loop of evolution by n+1/2 time step
  !! ***************************************************************************


  !! ***************************************************************************
  !! Exchange guardcell data for face-centered magnetic fields that will be
  !! used in the data reconstruction-evolution calculation.
  if ((.not. hy_forceHydroLimit) .and. hy_order > 1) then

     if (hy_fluxCorrect .and. hy_killdivb) then
        !! Call the electric field correction routine
        call Timers_start("conserveField")
        call Grid_conserveField()
        call Timers_stop("conserveField")
     endif

     !! Calculate divergence-free magnetic fields at interfaces
     do i=1,blockCount          !LOOP 2
        blockID = blockList(i)
        call Grid_getDeltas(blockID,del)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        halfTimeAdvance = .true.
        !! Evolve face-centered magnetic fields n+1/2 time step
        call hy_uhd_staggeredDivb(blockID,0.5*dt,del,blkLimits,blkLimitsGC,halfTimeAdvance)
     enddo

     !! Fill guardcells for only the face-centered field variable that are
     !! evolved to n+1/2 time step for date reconstruction
     gcMask  =.false.
#ifdef MAGI_FACE_VAR
     gcMask(MAGI_FACE_VAR:NFACE_VARS*NDIM:NFACE_VARS) = .true.
#endif
#ifdef MAG_FACE_VAR
     gcMask(MAG_FACE_VAR:NFACE_VARS*NDIM:NFACE_VARS) = .true. !DEV-DL: this could be false.
#endif

#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stamp('gcNeed(MAGI_FACE_VAR,MAG_FACE_VAR) - FACES', '[hy_uhd_unsplit]')
     end if
#endif

     !! Guardcell filling for facevars
     call Grid_fillGuardCells(FACES,ALLDIR,&
          maskSize=(hy_gcMaskSize-NUNK_VARS),mask=gcMask,makeMaskConsistent=.false.,&
          doLogMask=.NOT.gcMaskLogged)
  endif



  !! ***************************************************************************
  !! Second part of advancement to get Riemann states                          *
  !! ***************************************************************************
  !! Loop over the blocks
  do i=1,blockCount             !LOOP 3

     blockID = blockList(i)
     call Grid_getDeltas(blockID,del)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

     call hy_memGetBlkPtr(blockID,scrchFaceXPtr,SCRATCH_FACEX)
     if (NDIM > 1) call hy_memGetBlkPtr(blockID,scrchFaceYPtr,SCRATCH_FACEY)
     if (NDIM > 2) call hy_memGetBlkPtr(blockID,scrchFaceZPtr,SCRATCH_FACEZ)
     
#ifndef FIXEDBLOCKSIZE
     allocate(flx(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(fly(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(flz(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(    gravX(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(    gravY(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(    gravZ(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(  faceAreas(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif

     !! ************************************************************************
     !! Get gravity
     gravX = 0.
     gravY = 0.
     gravZ = 0.
     if (hy_useGravity) then
        ! We need to call gravity again here for cases (AMR or UG) WITH scratch vars when
        ! hy_useGravHalfUpdate=.false.
        ! This gravity is for the final update in the corrector step.
        call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
        gravX = gravX/hy_gref
        gravY = gravY/hy_gref
        gravZ = gravZ/hy_gref
     endif


     if (.not. hy_forceHydroLimit .and. hy_updateHydroFluxes .and. hy_order > 1) then
     !if (hy_updateHydroFluxes .and. hy_order > 1) then
        !! *********************************************************************
        !! Calculate Riemann (interface) states with divergence-free fields at n+1/2
        !! Note: gravX(:,:,:) - gravity at n
        call Timers_start("getRiemannState")
        call hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,dt,del, &
                                    gravX(:,:,:),gravY(:,:,:),gravZ(:,:,:),&
                                    scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,&
                                    normalFieldUpdate=.true.)
        call Timers_stop("getRiemannState")
     endif !end of if (hy_updateHydroFluxes) then

#if (NSPECIES+NMASS_SCALARS) > 0
     if (hy_fullSpecMsFluxHandling) then
        call hy_memGetBlkPtr(blockID,hy_spcR,CENTER)
        call hy_memGetBlkPtr(blockID,hy_spcL,SCRATCH_CTR)
     end if
#endif

#endif
! end if of #if NDIM > 1
#endif 
! end if of #if NFACE_VARS > 0
     !!************************************************************************************!!
     !!*** End of multidimensional MHD reconstruction (Riemann state) calculations ********!!
     !!************************************************************************************!!



     call hy_memGetBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)
     !! ************************************************************************
     !! Calculate high order Godunov fluxes
     !! Initialize arrays with zero
     flx = 0.
     fly = 0.
     flz = 0.
     call Timers_start("getFaceFlux")
     call hy_uhd_getFaceFlux(blockID,blkLimits,blkLimitsGC,datasize,del,flx,fly,flz,&
                             scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,scrch_Ptr,&
                             hy_SpcR,hy_SpcL,&
                             lastCall=.true.)
     call Timers_stop("getFaceFlux")

     if (hy_fullRiemannStateArrays) then
        call hy_memReleaseBlkPtr(blockID,scrchFaceXPtr,SCRATCH_FACEX)
        if (NDIM > 1) call hy_memReleaseBlkPtr(blockID,scrchFaceYPtr,SCRATCH_FACEY)
        if (NDIM > 2) call hy_memReleaseBlkPtr(blockID,scrchFaceZPtr,SCRATCH_FACEZ)
     else
        deallocate(scrchFaceXPtr)
        deallocate(scrchFaceYPtr)
        deallocate(scrchFaceZPtr)
     end if
      
#if (NSPECIES+NMASS_SCALARS) > 0
     if (hy_fullSpecMsFluxHandling) then
        if (hy_fullRiemannStateArrays .AND. NFACE_VARS > 0 .AND.  NDIM > 1) then
           call hy_memReleaseBlkPtr(blockID,hy_spcR,CENTER)
           call hy_memReleaseBlkPtr(blockID,hy_spcL,SCRATCH_CTR)
        else
           deallocate(hy_SpcR)
           deallocate(hy_SpcL)
        end if
     end if
#endif

     updateMode = UPDATE_NONE
#ifndef FLASH_GRID_UG
     if ((.not. hy_fullRiemannStateArrays) .OR. &
          (hy_fullSpecMsFluxHandling .AND. hy_numXN > 0) .OR. &
          .not. blockNeedsFluxCorrect(blockID)) then
#endif
        !! ************************************************************************
        !! Unsplit update for conservative variables from n to n+1 time step
        if (blockNeedsFluxCorrect(blockID)) then
           if (.NOT. hy_fullRiemannStateArrays) then
              updateMode = UPDATE_INTERIOR
           else
              updateMode = UPDATE_SPECMS_INTERIOR          
           end if
        else
           updateMode = UPDATE_ALL
        end if
#ifndef FLASH_GRID_UG
     end if
#endif
     if (updateMode == UPDATE_ALL) then ! done with updating all cells of this block after this loop iteration?
#if NFACE_VARS > 0
        !! ************************************************************************
        !! Update face-centered magnetic fields from n to n+1 time step
        if ((.not. hy_forceHydroLimit) .and. hy_killdivb .and. hy_order > 1) then
           halfTimeAdvance = .false.
           call hy_uhd_getElectricFields(blockID,blkLimits,blkLimitsGC,del,flx,fly,flz)
! <- ychen 09-2014
           if (dr_simTime.ge.sim(nozzle)%tOn .and. &
               dr_simTime.lt.(sim(nozzle)%tOn+sim(nozzle)%duration)) then
              call hy_uhd_electricNozzle(blockID,blkLimits,blkLimitsGC)
           endif
! ychen ->
           call hy_uhd_staggeredDivb(blockID,dt,del,blkLimits,blkLimitsGC,halfTimeAdvance)
        endif ! End of if ((.not. hy_forceHydroLimit) .and. hy_killdivb .and. hy_order > 1) then
#endif

        !! ************************************************************************
        !! Unsplit update for conservative variables from n to n+1 time step
        call Timers_start("unsplitUpdate")
        call hy_uhd_unsplitUpdate(blockID,updateMode,dt,del,datasize,blkLimits,&
                                  blkLimitsGC,flx,fly,flz,gravX,gravY,gravZ,&
                                  scrch_Ptr)
        call Timers_stop("unsplitUpdate")

#ifndef GRAVITY /* if gravity is included we delay energy fix until we update gravity at n+1 state */
        !! Correct energy if necessary
        call hy_uhd_energyFix(blockID,blkLimits,dt,dtOld,del,hy_unsplitEosMode)

#ifdef FLASH_UHD_3T
        call Timers_start("unsplitUpdate 3T")
        call hy_uhd_unsplitUpdateMultiTemp&
             (blockID,updateMode,blkLimits, dataSize, dt, del, flx, fly, flz, scrch_Ptr)
        call Timers_stop("unsplitUpdate 3T")
#endif

        if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
           !! Convert unit
           call hy_uhd_unitConvert(blockID,BWDCONVERT)
        endif

!#ifndef FLASH_EOS_GAMMA
        !! Call to Eos
        call Eos_wrapped(hy_eosModeAfter, blkLimits, blockID)
!#endif
#endif /* ifndef GRAVITY */
     else if (updateMode==UPDATE_INTERIOR) then
        !! Update interior cells now, defer update of one layer of cells until later.

        call Timers_start("unsplitUpdate")
        call hy_uhd_unsplitUpdate(blockID,updateMode,dt,del,datasize,blkLimits,&
                                  blkLimitsGC,flx,fly,flz,gravX,gravY,gravZ,&
                                  scrch_Ptr)
        call Timers_stop("unsplitUpdate")

#ifdef FLASH_UHD_3T
        call Timers_start("unsplitUpdate 3T")
        call hy_uhd_unsplitUpdateMultiTemp&
                (blockID,updateMode,blkLimits, dataSize, dt, del, flx, fly, flz, scrch_Ptr)
        call Timers_stop("unsplitUpdate 3T")
#endif

     else if (updateMode==UPDATE_SPECMS_INTERIOR) then
        !! if (nearly) everything is left to do after Flux correction.
        !! Flux conservation calls on AMR:
        !! Correct fluxes at each block boundary where coarse and fine
        !! blocks are neighboring each other.

#if NFACE_VARS > 0
        !! ************************************************************************
        !! Update face-centered magnetic fields from n to n+1 time step
        if ((.not. hy_forceHydroLimit) .and. hy_killdivb .and. hy_order > 1) then
           halfTimeAdvance = .false.
           call hy_uhd_getElectricFields(blockID,blkLimits,blkLimitsGC,del,flx,fly,flz)
! <- ychen 09-2014
  !write(*,*) "*************** electricNozzle1 ***************"
           if (dr_simTime.ge.sim(nozzle)%tOn .and.&
               dr_simTime.lt.(sim(nozzle)%tOn+sim(nozzle)%duration)) then
              call hy_uhd_electricNozzle(blockID,blkLimits,blkLimitsGC)
           endif
  !write(*,*) "*************** electricNozzle2 ***************"
! ychen ->
           !write(*,*) blockID, 'unsplit'
        endif ! End of if ((.not. hy_forceHydroLimit) .and. hy_killdivb .and. hy_order > 1) then
#endif
        call Timers_start("unsplitUpdate")
        call hy_uhd_unsplitUpdate(blockID,updateMode,dt,del,datasize,blkLimits,&
                                  blkLimitsGC,flx,fly,flz,gravX,gravY,gravZ,&
                                  scrch_Ptr)
        call Timers_stop("unsplitUpdate")

     else if (updateMode==UPDATE_NONE) then
        !! if (nearly) everything is left to do after Flux correction.
        !! Flux conservation calls on AMR:
        !! Correct fluxes at each block boundary where coarse and fine
        !! blocks are neighboring each other.

#if NFACE_VARS > 0
        !! ************************************************************************
        !! Update face-centered magnetic fields from n to n+1 time step
        if ((.not. hy_forceHydroLimit) .and. hy_killdivb .and. hy_order > 1) then
           halfTimeAdvance = .false.
           call hy_uhd_getElectricFields(blockID,blkLimits,blkLimitsGC,del,flx,fly,flz)
        endif ! End of if ((.not. hy_forceHydroLimit) .and. hy_killdivb .and. hy_order > 1) then
#endif
     end if

     if (blockMustStoreFluxes(blockID)) then
        if (hy_geometry /= CARTESIAN) then
           ! we are using consv_fluxes and need to divide by face areas
           call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, &
                                (/1,1,1/), faceAreas, datasize)
           call Grid_putFluxData(blockID,IAXIS,flx,datasize,hy_fluxCorVars,faceAreas)

           if (NDIM > 1) then
              faceAreas = 0.
              call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
                                   (/1,1,1/), faceAreas, datasize)
              call Grid_putFluxData(blockID,JAXIS,fly,datasize,hy_fluxCorVars,faceAreas)

              if (NDIM > 2) then
                 faceAreas = 0.
                 call Grid_getBlkData(blockID,CELL_FACEAREA,KLO_FACE, EXTERIOR, &
                                      (/1,1,1/), faceAreas, datasize)
                 call Grid_putFluxData(blockID,KAXIS,flz,datasize,hy_fluxCorVars,faceAreas)
              endif
           endif
        else ! Cartesian geometry
           call Grid_putFluxData(blockID,IAXIS,flx,datasize)
           if (NDIM > 1) then
              call Grid_putFluxData(blockID,JAXIS,fly,datasize)
              if (NDIM > 2) then
                 call Grid_putFluxData(blockID,KAXIS,flz,datasize)
              endif
           endif
        endif

     endif

     
#ifndef FIXEDBLOCKSIZE
     deallocate(flx)
     deallocate(fly)
     deallocate(flz)
     deallocate(gravX)
     deallocate(gravY)
     deallocate(gravZ)
     deallocate(faceAreas)
#endif

     call hy_memReleaseBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)
  end do
  !! End of leaf block do-loop before flux conserve call

#if (NSPECIES+NMASS_SCALARS) > 0
  if (hy_fullSpecMsFluxHandling .AND. hy_fullRiemannStateArrays &
       .AND. NFACE_VARS > 0 .AND. NDIM > 1) then
     call hy_memDeallocScratch(CENTER,      arrayRank=5)
     call hy_memDeallocScratch(SCRATCH_CTR, arrayRank=5)
  end if
#endif


  !! ***************************************************************************
  !! Third part of advancement                                                 *
  !! ***************************************************************************
  !! Do this part only if refining and flux correcting
  if (hy_fluxCorrect) then

     !! ************************************************************************
     !! Conservation of Fluxes at each block boundary
     call Timers_start("conserveFluxes")
     call Grid_conserveFluxes(ALLDIR,level)
     call Timers_stop("conserveFluxes")

#if NFACE_VARS > 0
     !! ************************************************************************
     !! Correct electric fields at fine-coarse boundaries
     if ((.not. hy_forceHydroLimit) .and. (hy_killdivb)) then
        !! Call the electric field correction routine
        call Timers_start("conserveField")
        call Grid_conserveField()
        call Timers_stop("conserveField")
     endif
#endif


     !! ************************************************************************
     !! Perform updates for every blocks
     do i = 1,blockCount        !LOOP 4

        blockID = blockList(i)
        if (.NOT. blockNeedsFluxCorrect(blockID)) CYCLE ! Save a lot of work!

        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        call Grid_getDeltas(blockID,del)

        datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
#ifndef FIXEDBLOCKSIZE
        allocate(flx(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(fly(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(flz(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(    gravX(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(    gravY(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(    gravZ(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(  faceAreas(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif

        !! *********************************************************************
        !! Get gravity
        gravX = 0.
        gravY = 0.
        gravZ = 0.
        if (hy_useGravity) then
           call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
           gravX = gravX/hy_gref
           gravY = gravY/hy_gref
           gravZ = gravZ/hy_gref
        endif

        !! Get modified conserved flux values at block interfaces:
        !! This is important especially at the block interface at 
        !! different refinement levels of neighboring blocks
        flx = 0.
        fly = 0.
        flz = 0.

        if (hy_fullRiemannStateArrays) then
           call hy_memGetBlkPtr(blockID,scrchFaceXPtr,SCRATCH_FACEX)
           if (NDIM > 1) call hy_memGetBlkPtr(blockID,scrchFaceYPtr,SCRATCH_FACEY)
           if (NDIM > 2) call hy_memGetBlkPtr(blockID,scrchFaceZPtr,SCRATCH_FACEZ)
        endif

        if (hy_fullRiemannStateArrays) then
           flx(HY_DENS_FLUX:HY_END_FLUX,&
                blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))&
              =scrchFaceXPtr(HY_P01_FACEXPTR_VAR:HY_P01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,&
                         blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
                         blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                         blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
#if NDIM > 1
           fly(HY_DENS_FLUX:HY_END_FLUX,&
                blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,&
                blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))&
              =scrchFaceYPtr(HY_P01_FACEYPTR_VAR:HY_P01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,&
                         blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                         blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,&
                         blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
#if NDIM == 3
           flz(HY_DENS_FLUX:HY_END_FLUX,&
                blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)&
              =scrchFaceZPtr(HY_P01_FACEZPTR_VAR:HY_P01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,&
                         blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                         blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                         blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)
#endif
#endif
        endif

        if (hy_geometry /= CARTESIAN) then
           call Grid_getBlkData(blockID,CELL_FACEAREA,ILO_FACE, EXTERIOR, &
                                (/1,1,1/), faceAreas, datasize)
           call Grid_getFluxData(blockID,IAXIS,flx,datasize,hy_fluxCorVars,faceAreas)
        else
           call Grid_getFluxData(blockID,IAXIS,flx,datasize)
        endif

#if NDIM > 1
        if (hy_geometry /= CARTESIAN) then
           call Grid_getBlkData(blockID,CELL_FACEAREA,JLO_FACE, EXTERIOR, &
                                (/1,1,1/), faceAreas, datasize)
           call Grid_getFluxData(blockID,JAXIS,fly,datasize,hy_fluxCorVars,faceAreas)
        else
           call Grid_getFluxData(blockID,JAXIS,fly,datasize)
        endif
#if NDIM > 2
        if (hy_geometry /= CARTESIAN) then
           call Grid_getBlkData(blockID,CELL_FACEAREA,KLO_FACE, EXTERIOR, &
                                (/1,1,1/), faceAreas, datasize)
           call Grid_getFluxData(blockID,KAXIS,flz,datasize,hy_fluxCorVars,faceAreas)
        else
           call Grid_getFluxData(blockID,KAXIS,flz,datasize)
        endif
#endif
#endif

        if (hy_fullRiemannStateArrays) then
           call hy_memReleaseBlkPtr(blockID,scrchFaceXPtr,SCRATCH_FACEX)
           if (NDIM > 1) call hy_memReleaseBlkPtr(blockID,scrchFaceYPtr,SCRATCH_FACEY)
           if (NDIM > 2) call hy_memReleaseBlkPtr(blockID,scrchFaceZPtr,SCRATCH_FACEZ)
        endif

        call hy_memGetBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)
        !! *********************************************************************
        !! Unsplit update for conservative variables from n to n+1 time step
        if (.not. hy_fullRiemannStateArrays) then
           updateMode=UPDATE_BOUND
        else if (hy_fullSpecMsFluxHandling .AND. hy_numXN > 0) then
           updateMode=UPDATE_ALL_SPECMS_BOUND
        else
           updateMode=UPDATE_ALL
        endif

#if NFACE_VARS > 0
        !! *********************************************************************
        !! Update face-centered magnetic fields from n to n+1 time step
        if (.not. hy_forceHydroLimit .and. hy_killdivb) then
           halfTimeAdvance = .false.
           call hy_uhd_staggeredDivb(blockID,dt,del,blkLimits,blkLimitsGC,halfTimeAdvance)
        endif
#endif

        !! *********************************************************************
        !! Unsplit update for conservative variables from n to n+1 time step
        call Timers_start("unsplitUpdate")
        call hy_uhd_unsplitUpdate(blockID,updateMode,dt,del,datasize,blkLimits,&
                                  blkLimitsGC,flx,fly,flz,gravX,gravY,gravZ,&
                                  scrch_Ptr)
        call Timers_stop("unsplitUpdate")


#ifndef GRAVITY /* if gravity is included we delay energy fix until we update gravity at n+1 state */
        !! *********************************************************************
        !! Correct energy if necessary
        call hy_uhd_energyFix(blockID,blkLimits,dt,dtOld,del,hy_unsplitEosMode)

#ifdef FLASH_UHD_3T
        call Timers_start("unsplitUpdate 3T")
        call hy_uhd_unsplitUpdateMultiTemp&
             (blockID,updateMode,blkLimits, dataSize, dt, del, flx, fly, flz, scrch_Ptr)
        call Timers_stop("unsplitUpdate 3T")
#endif

        !! *********************************************************************
        !! Convert units if necessary
        if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
           !! Convert units
           call hy_uhd_unitConvert(blockID,BWDCONVERT)
        endif


!#ifndef FLASH_EOS_GAMMA
        !! Call to Eos
        call Eos_wrapped(hy_eosModeAfter, blkLimits, blockID)
!#endif
#endif /* if gravity is included, we calculate gravity for n+1 in the below */

        call hy_memReleaseBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)

#ifndef FIXEDBLOCKSIZE
        deallocate(flx)
        deallocate(fly)
        deallocate(flz)
        deallocate(gravX)
        deallocate(gravY)
        deallocate(gravZ)
        deallocate(faceAreas)
#endif
        
     end do !! End of the loop over blocks
  end if !! End of the flux conservation routine

  call Timers_start("tail")

  if (hy_fullRiemannStateArrays) then
     call hy_memDeallocScratch(SCRATCH_FACEX)
     if (NDIM > 1) call hy_memDeallocScratch(SCRATCH_FACEY)
     if (NDIM > 2) call hy_memDeallocScratch(SCRATCH_FACEZ)
  endif


#ifdef GRAVITY /* Perform this only when gravity is used */
  !! ***************************************************************************
  !! Fourth part of advancement to compute gravity at n+1 state                *
  !! ***************************************************************************

#ifdef GPOT_VAR
  if (hy_useGravity) then
     ! The following call invokes Gravity_potentialListOfBlocks and related stuff,
     ! to prepare for retrieving updated accelerations below.
     call hy_uhd_prepareNewGravityAccel(blockCount,blockList,gcMaskLogged)
  endif
#endif


  !! Proceed to couple the updated gravitational accelerations 
  !! to energy and momenta on each block (see hy_uhd_addGravityUnsplit)
  do i = 1,blockCount           !LOOP 5

     blockID = blockList(i)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getDeltas(blockID,del)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
#ifndef FIXEDBLOCKSIZE
     allocate(gravX(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravY(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravZ(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif

     !! *********************************************************************
     !! Get and add gravity to hydro variables (momenta and energy)         *
     !! *********************************************************************
     gravX = 0.
     gravY = 0.
     gravZ = 0.
     if (hy_useGravity) then

        call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ,&
             lastCall=.TRUE.)
        gravX = gravX/hy_gref
        gravY = gravY/hy_gref
        gravZ = gravZ/hy_gref

        call hy_uhd_addGravityUnsplit(blockID,blkLimits,dataSize,dt,&
             gravX(:,:,:),gravY(:,:,:),gravZ(:,:,:))
     endif

     !! *********************************************************************
     !! Correct energy if necessary
     call hy_uhd_energyFix(blockID,blkLimits,dt,dtOld,del,hy_unsplitEosMode)

     !! *********************************************************************
     !! Convert units if necessary
     if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
        !! Convert units
        call hy_uhd_unitConvert(blockID,BWDCONVERT)
     endif

     !! *********************************************************************
     !! Call to Eos
     call Eos_wrapped(hy_eosModeAfter, blkLimits, blockID)

#ifndef FIXEDBLOCKSIZE
     deallocate(gravX)
     deallocate(gravY)
     deallocate(gravZ)
#endif

  end do !! End of the loop over blocks
#endif /* End of n+1 gravity coupling */

  ! Do 3T update...
#ifdef FLASH_UHD_3T
  call hy_uhd_multiTempAfter(blockCount, blockList, dt)
#endif



#if NFACE_VARS > 0
  !! ***********************************************************************
  !! Note: Need guardcell fills to compute current densities and divv.
  !!       Use this call only if needed as this is expensive after all!
#ifdef HY_NEED_EXTRA_GCFILL


#ifdef FLASH_GRID_PARAMESH
  gcMask  =.FALSE.
  gcMask(MAGX_VAR)=.true.
  gcMask(MAGY_VAR)=.true.
  gcMask(MAGZ_VAR)=.true.
  gcMask(VELX_VAR)=.true.
  gcMask(VELY_VAR)=.true.
  gcMask(VELZ_VAR)=.true.
#ifdef MAG_FACE_VAR
  gcMask(NUNK_VARS+MAG_FACE_VAR:NUNK_VARS+NFACE_VARS*NDIM:NFACE_VARS) = .true.
#endif

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     call Logfile_stampVarMask(gcMask, .FALSE., '[hy_uhd_unsplit]', 'gcWant')
  end if
#endif
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,minLayers=1,doEos=.false.,&
       maskSize=hy_gcMaskSize,mask=gcMask,makeMaskConsistent=.true.,&
       doLogMask=.NOT.gcMaskLogged)
#endif

#ifdef FLASH_GRID_UG
  !! Note that, on UG, guardcell masking is not applicable
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR)
#endif


  !! DEV: DL-there should be cylindrical implementations in differencings here!
  do i=1,blockCount             !LOOP 6

     blockID = blockList(i)
     call Grid_getBlkPtr(blockID,U,CENTER)
     call Grid_getDeltas(blockID,del)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     do iz=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do iy=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do ix=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

              !! current density jx
#ifdef CURX_VAR
#if NDIM > 1
              U(CURX_VAR,ix,iy,iz)=&
                   0.5*(U(MAGZ_VAR,ix,iy+1,iz)-U(MAGZ_VAR,ix,iy-1,iz))/del(DIR_Y)
#if NDIM > 2
              U(CURX_VAR,ix,iy,iz)=U(CURX_VAR,ix,iy,iz)&
                   -0.5*(U(MAGY_VAR,ix,iy,iz+1)-U(MAGY_VAR,ix,iy,iz-1))/del(DIR_Z)
#endif
#endif
#endif

              !! current density jy
#ifdef CURY_VAR
              U(CURY_VAR,ix,iy,iz)=&
                   -0.5*(U(MAGZ_VAR,ix+1,iy,iz)-U(MAGZ_VAR,ix-1,iy,iz))/del(DIR_X)
#if NDIM > 2
              U(CURY_VAR,ix,iy,iz)=U(CURY_VAR,ix,iy,iz)&
                   +0.5*(U(MAGX_VAR,ix,iy,iz+1)-U(MAGX_VAR,ix,iy,iz-1))/del(DIR_Z)
#endif
#endif

              !! current density jz
#ifdef CURZ_VAR
              U(CURZ_VAR,ix,iy,iz)=&
                   0.5*(U(MAGY_VAR,ix+1,iy,iz)-U(MAGY_VAR,ix-1,iy,iz))/del(DIR_X)

#if NDIM > 1
              U(CURZ_VAR,ix,iy,iz)=U(CURZ_VAR,ix,iy,iz)&
                   -0.5*(U(MAGX_VAR,ix,iy+1,iz)-U(MAGX_VAR,ix,iy-1,iz))/del(DIR_Y)
#endif
#endif

#ifdef DIVV_VAR

              U(DIVV_VAR,ix,iy,iz) = (U(VELX_VAR,ix+1,iy,iz)-U(VELX_VAR,ix-1,iy,iz))/del(DIR_X)
              if (NDIM > 1) then
                 U(DIVV_VAR,ix,iy,iz) = U(DIVV_VAR,ix,iy,iz) + &
                                   +(U(VELY_VAR,ix,iy+1,iz)-U(VELY_VAR,ix,iy-1,iz))/del(DIR_Y)
              if (NDIM == 3) then
                 U(DIVV_VAR,ix,iy,iz) = U(DIVV_VAR,ix,iy,iz) + &
                                   +(U(VELZ_VAR,ix,iy,iz+1)-U(VELZ_VAR,ix,iy,iz-1))/del(DIR_Z)
              endif
              endif
#endif
           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blockID,U,CENTER)
  enddo


#endif 
#endif


  ! ***************************************
  ! *                                     *
  ! *     ADD BIERMANN BATTERY SOURCE     *
  ! *                                     *
  ! ***************************************
  !
  ! This subroutine call modifies the magnetic field directly using an
  ! external Biermann Battery source term. It will only preserve
  ! div(B)=0 in xy and rz geometries ...
  if ( (hy_useBiermann .or. hy_useBiermann1T) .and. hy_biermannSource ) then
#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stamp('gcNeed - CENTER unmasked', '[hy_uhd_unsplit]')
     end if
#endif
     call Grid_fillGuardCells(CENTER, ALLDIR)
     call hy_uhd_biermannSource(blockCount, blockList, dt)
  end if

  call hy_memDeallocScratch()

  call Driver_getSimTime(hy_simTime, hy_simGeneration)


#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     gcMaskLogged = .TRUE.
  end if
#endif

  call Timers_stop("tail")

End Subroutine hy_uhd_unsplit
