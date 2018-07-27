!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_unsplitUpdate
!!
!! NAME
!!
!!  hy_uhd_unsplitUpdate
!!
!! SYNOPSIS
!!
!!  call hy_uhd_unsplitUpdate( integer(IN) :: blockID,
!!                        integer(IN) :: rangeSwitch,
!!                        real(IN)    :: dt,
!!                        real(IN)    :: del(MDM),
!!                        integer(IN) :: dataSize(MDIM),
!!                        integer(IN) :: blkLimits(2,MDIM),
!!                        integer(IN) :: blkLimitsGC(2,MDIM),
!!                        real(IN)    :: xflux(:,:,:,:), 
!!                        real(IN)    :: yflux(:,:,:,:), 
!!                        real(IN)    :: zflux(:,:,:,:),
!!                        real(IN)    :: gravX(:,:,:),
!!                        real(IN)    :: gravY(:,:,:),
!!                        real(IN)    :: gravZ(:,:,:),
!!                        real,POINTER, dimension(:,:,:,:) :: scrch_Ptr)
!!
!! ARGUMENTS
!!
!!   blockID      - current block ID
!!   rangeSwitch  - switch for selective updates on AMR (only valid for hydro and not for MHD).
!!                  One of the following values (defined in UHD.h)
!!                    UPDATE_ALL,              update all variables in all cells.
!!                    UPDATE_INTERIOR,         update all variables in interior cells.
!!                    UPDATE_BOUND,            update all variables in boundary cells.
!!                    UPDATE_SPECMS_INTERIOR,  update species and mass scalar variables in interior cells;
!!                                             values of those variables in those cells
!!                                             will be left in multiplied-by-density form.
!!                    UPDATE_ALL_SPECMSBOUND,  update non-species / non-mass-scalar variables in all
!!                                             cells; update species and mass scalar variables in
!!                                             boundary cells; and finish updating of species and mass
!!                                             scalar variables in interior cells,
!!                                             by dividing those values by the updated density.
!!                  Boundary cells here are cells that have at least one face at a block boundary;
!!                  interior cells are the remaining cells of a block.
!!                  Note that neither boundary cells nor interior cells, in the sense used here,
!!                  are guard cells! Guard cells should not be modified by this routine.
!!   dt           - timestep
!!   del          - coordinate spacings in {X,Y,Z} directions
!!   dataSize     - size of the current block
!!   blkLimits    - an array that holds the lower and upper indices of the section
!!                  of block without the guard cells
!!   blkLimitsGC    - an array that holds the lower and upper indices of the section
!!                  of block with the guard cells
!!   xflux,yflux,zflux - cell face centered fluxes at each {=X,Y,Z} direction
!!   gravX,gravY,gravZ - gravitational acceleration components in X,Y,Z directions
!!   scrch_Ptr      - pointer to scratch space; used both for receiving some data, and
!!                    for storing some data for later processing in the FLASH_UHD_3T case.
!!
!! DESCRIPTION
!!
!!   This routine updates the cell-centered conservative variables and intermediate
!!   internal energy to the next time step using fluxes for all spatial directions (unsplit scheme).
!!
!! NOTES
!!
!!   The ENER_VAR and EINT_VAR components of UNK are now always kept in specific
!!   (i.e., energy per mass) form, on entry as well as on return from this routine.
!!   This is changed from the behavior in FLASH4.2.2 and earlier.
!!***

!!REORDER(4): U, Uold, SpOld, scrch_Ptr, [xyz]flux

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "UHD.h"

  Subroutine hy_uhd_unsplitUpdate(blockID,rangeSwitch,dt,del,dataSize,blkLimits,&
                                  blkLimitsGC,xflux,yflux,zflux,gravX,gravY,gravZ,&
                                  scrch_Ptr)

    use Hydro_data,           ONLY : hy_smalldens,hy_order,hy_irenorm,hy_numXN, &
                                     hy_meshMe, &
                                     hy_geometry, hy_gcMaskSize, &
                                     hy_gcMask, hy_eswitch, hy_useAuxEintEqn, &
                                     hy_fullSpecMsFluxHandling, &
                                     hy_conserveAngMom, hy_smallE, hy_irenorm
    use Hydro_data,           ONLY : fP => hy_fPresInMomFlux
#ifdef FLASH_USM_MHD
  use Hydro_data,             ONLY : hy_mref, hy_hallVelocity, &
                                     hy_useMagneticResistivity, hy_conserveAngField
  use MagneticResistivity_interface, &
                              ONLY : MagneticResistivity
  use hy_uhd_interface,       ONLY : hy_uhd_addOhmicHeating 
#endif
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) 
    use Hydro_data,           ONLY : hy_forceHydroLimit
#endif
    use hy_uhd_interface,     ONLY : hy_uhd_updateSpeciesMassScalar
#ifdef FLASH_UHD_3T
#ifdef FLASH_USM_MHD
  use hy_uhd_interface,       ONLY : hy_uhd_getCurrents
#endif
#endif
    use Grid_interface,       ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
                                     Grid_getCellCoords, Grid_getBlkData,&
                                     Grid_renormAbundance, Grid_limitAbundance
#ifdef FLASH_UHD_3T
#ifdef FLASH_USM_MHD
    use Eos_interface,        ONLY : Eos_getAbarZbar
#endif
#endif

    implicit none

    !! ---- Arguments ---------------------------------
    integer,intent(IN) :: blockID, rangeSwitch
    real, intent(IN)   :: dt
    real, intent(IN)   :: del(MDIM)
    integer,dimension(MDIM),intent(IN) :: dataSize
    integer,intent(IN) :: blkLimits(LOW:HIGH,MDIM)
    integer,intent(IN) :: blkLimitsGC(LOW:HIGH,MDIM)
#ifdef FIXEDBLOCKSIZE
    real, intent(in) :: xflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
    real, intent(in) :: yflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
    real, intent(in) :: zflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
    real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: gravX,gravY,gravZ

#else
    real, intent(in) :: xflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS))
    real, intent(in) :: yflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS))  
    real, intent(in) :: zflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS))
    real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)),& 
         intent(IN) :: gravX,gravY,gravZ
#endif
    !!---------------------------------------------------

    integer :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
    real    :: dx, dy, dz, dx_sph
    real    :: dxv, dyv
    real, dimension(HY_VARINUM) :: U0
    real, dimension(NFLUXES)    :: FL,FR,GL,GR,HL,HR
    real, pointer, dimension(:,:,:,:) :: U, scrch_Ptr
    real    :: densNew, densNph ! densNph: "dens at n+1/2"
    real    :: IntEner,tempPres

#if (NSPECIES+NMASS_SCALARS) > 0
#ifdef FIXEDBLOCKSIZE
    real, dimension(hy_numXN,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: SpOld
    real, dimension(6,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: Uold
#else
    real, dimension(hy_numXN,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)) :: SpOld
    real, dimension(6,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)) :: Uold
#endif
#else
  !This is just here so I can have the same omp parallel directive whether
  !(NSPECIES+NMASS_SCALARS) > 0 or not.  SpOld and UOld are not used.
  real, dimension(1,1,1,1) :: SpOld
  real, dimension(1,1,1,1) :: UOld
#endif

    integer :: iSize, jSize, kSize
#ifdef FIXEDBLOCKSIZE
    real, dimension(GRID_IHI_GC) :: xCenter, xLeft, xRight
#else
    real, dimension(dataSize(IAXIS)) :: xCenter, xLeft, xRight
#endif
#if NDIM == 3
#  ifdef FIXEDBLOCKSIZE
    real, dimension(GRID_JLO:GRID_JHI) :: yCenter
#  else
    real, dimension(blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)) :: yCenter
#  endif
#else
    real, dimension(0) :: yCenter
#endif

    logical :: gcMask(hy_gcMaskSize)
    real, dimension(HY_VARINUM) :: Sgeo, Sphys
    real, dimension(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
                    blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),  &
                    blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS) ) :: faceAreas
#if NDIM > 1
    real, dimension(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)  ,&
                    blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,&
                    blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS) ) :: faceAreasY
#endif
    real, dimension(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                    blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),  &
                    blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS) ) :: cellVolumes

#if defined(FLASH_USM_MHD) && defined(FLASH_UHD_3T) 
    real, dimension(3, dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)) :: Jp, Jm
#else
    real    :: Jp, Jm
#endif
    integer :: ky,kz,iskip
    real    :: temp, hdt
    real    :: presStar, densStar, pmomStar, tmomStar, xmomStar
    real    :: pmagStar, xmagStar, zmagStar
    real    :: xvel0, phiVel0, thtVel0
    integer :: VEL_PHI, MOM_PHI, MOM_PHI_FLUX, MAG_PHI,  MAG_PHI_FLUX
    integer :: VEL_ZI, MOM_ZI, MOM_ZI_FLUX, MAG_ZI,  MAG_ZI_FLUX
    integer :: VEL_THT, MOM_THT, MOM_THT_FLUX, velPhiVar
    real    :: leftFac, rghtFac, dPdr, rvol, alpha, cs, eta
    real    :: ekin, eint, newEint
    real    :: abar, zbar
    real    :: Qohm
#if (NSPECIES+NMASS_SCALARS) > 0
    integer :: isph,ispu
    real    :: sumSpecies
    real    :: newDens
    logical :: specialForInterior
#endif
  !! Resistive MHD 
#ifdef FIXEDBLOCKSIZE 
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, & 
                  GRID_JLO_GC:GRID_JHI_GC, & 
                  GRID_KLO_GC:GRID_KHI_GC) :: res_eta

  real, dimension(GRID_ILO_GC:GRID_IHI_GC, & 
                  GRID_JLO_GC:GRID_JHI_GC, & 
                  GRID_KLO_GC:GRID_KHI_GC) :: res_source
                  
#else 
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), & 
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), & 
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: res_eta

  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), & 
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), & 
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: res_source                
#endif 
  real    :: eta_loc, Jyp, Jym, dxBzm,dxBzp, inv_dVrm, inv_dVrp
  real, pointer,dimension(:)    :: speciesArr
#ifdef FLASH_USM_MHD
  res_eta = 0.0
#endif


    hdt = 0.5*dt

    iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
    jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
    kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

    if (rangeSwitch==UPDATE_SPECMS_INTERIOR) then
       call unsplitUpdateSpecMs
       return                   ! We are already done, RETURN now.
    end if

    Jp = 0.0 
    Jm = 0.0

    !! Set ranges for update
    imin  = blkLimits(LOW, IAXIS)
    imax  = blkLimits(HIGH,IAXIS)
    jmin  = 1
    jmax  = 1
    kmin  = 1
    kmax  = 1

    dx = del(DIR_X)
    dy = 1.
    dz = 1.
    if (NDIM >= 2) then
       jmin  = blkLimits(LOW, JAXIS)
       jmax  = blkLimits(HIGH,JAXIS)
       dy = del(DIR_Y)
       if (NDIM == 3) then
          kmin  = blkLimits(LOW, KAXIS)
          kmax  = blkLimits(HIGH,KAXIS)
          dz = del(DIR_Z)
       endif
    endif
    dxv = dx

    !! Set regions to update depending on update mode
    iskip = 1
    if (rangeSwitch==UPDATE_INTERIOR) then
       imin  = imin+1
       imax  = imax-1
       !iskip = 1
       if (NDIM >= 2) then
          jmin  = jmin+1
          jmax  = jmax-1
          if (NDIM == 3) then
             kmin  = kmin+1
             kmax  = kmax-1
          endif
       endif
    endif


    ! Get block pointers
    call Grid_getBlkPtr(blockID,U,CENTER)

    if (hy_geometry /= CARTESIAN) then
       faceAreas = 0.
       call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, &
            (/blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS)/), &
            faceAreas(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),  &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)), &
            (/isize+1, jsize, ksize/) )
#if NDIM > 1
       if (hy_geometry == SPHERICAL) then
          call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
            (/blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS)/), &
            faceAreasY(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,  &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)), &
            (/isize, jsize+1, ksize/) )
       end if
#endif

       cellVolumes = 0.
       call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR, &
            (/blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS)/), &
            cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)), &
            (/isize, jsize, ksize/) )
    endif

    
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
    if (hy_forceHydroLimit) then
       U(MAGX_VAR:MAGZ_VAR,:,:,:) = 0.
    endif
#endif


#if (NSPECIES+NMASS_SCALARS) > 0
    do ispu =  SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
       isph= ispu-NPROP_VARS
       SpOld(isph,:,:,:) = U(ispu,:,:,:)
    enddo
    Uold(1,  :,:,:) = U(DENS_VAR,:,:,:)
    Uold(2,  :,:,:) = U(PRES_VAR,:,:,:)
    Uold(3:5,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)
    Uold(6,  :,:,:) = U(GAME_VAR,:,:,:)
#endif

#ifdef DEBUG_HYDRO_POSITIVITY
    call    Grid_getCellCoords(IAXIS,blockID, CENTER,    .true.,xCenter, dataSize(IAXIS))
#else
    if (hy_geometry /= CARTESIAN) then
       call Grid_getCellCoords(IAXIS,blockID, CENTER,    .true.,xCenter, dataSize(IAXIS))
    end if
#endif
    if (hy_geometry /= CARTESIAN) then
       call Grid_getCellCoords(IAXIS,blockID, LEFT_EDGE, .true.,xLeft,   dataSize(IAXIS))
       call Grid_getCellCoords(IAXIS,blockID, RIGHT_EDGE,.true.,xRight,  dataSize(IAXIS))
       if (NDIM == 3 .AND. hy_geometry == SPHERICAL) then
          call Grid_getCellCoords(JAXIS,blockID, CENTER,.false.,yCenter, size(yCenter))
       end if
    endif

#ifdef FLASH_USM_MHD
#ifdef FLASH_UHD_3T
    !! STORE THE OLD FIELD FOR CURRENT CALCULATION FOR 3T UPDATE
    scrch_Ptr(HY_XN01_SCRATCHCTR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = U(MAGX_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
    scrch_Ptr(HY_XN02_SCRATCHCTR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = U(MAGY_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
    scrch_Ptr(HY_XN03_SCRATCHCTR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = U(MAGZ_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
#endif
#endif

    !! Define dimension dependent switches
    ky=0
    kz=0

    if (NDIM > 1) then
       ky=1
       if (NDIM > 2) then
          kz=1
       endif
    endif


    iskip = 1
    if (NDIM == 1 .and. rangeSwitch .eq. UPDATE_BOUND) iskip = imax-imin
    !! Loop to get magnetic resistivity source correction 
#ifdef FLASH_USM_MHD
    if (hy_geometry == CYLINDRICAL .and. hy_useMagneticResistivity) then
       do k=kmin,kmax
          do j=jmin-2,jmax+2
             if (NDIM >= 2) then
                iskip = 1
                if (rangeSwitch == UPDATE_BOUND .and. j > jmin .and. j < jmax) then
                   iskip = imax-imin
                   if (NDIM == 3) then
                      iskip = 1
                      if (k > kmin .and. k < kmax) then
                         iskip = imax-imin
                      endif
                   endif
                endif
             endif

             do i=imin-2,imax+2,iskip
                !! Get magnetic Resistivity
                call MagneticResistivity(U(:,i-1,j,k),res_eta(i-1,j,k))
                call MagneticResistivity(U(:,i,j,k),res_eta(i,j,k))
                call MagneticResistivity(U(:,i+1,j,k),res_eta(i+1,j,k))                  

                !! Normalize if needed
                res_eta(i-1,j,k)  = res_eta(i-1,j,k)/hy_mref
                res_eta(i  ,j,k)  = res_eta(i  ,j,k)/hy_mref
                res_eta(i+1,j,k)  = res_eta(i+1,j,k)/hy_mref

                !! We are adding the d (eta Jz)/dr as a source
                !! term so as to keep the induction of Bphi with the source
                !! formulation as it does not balance without....
                eta_loc = 0.5*(res_eta(i,j,k)+res_eta(i-1,j,k))
                inv_dVrm = xCenter(i)*abs(xCenter(i)) - xCenter(i-1)*abs(xCenter(i-1))
                inv_dVrm = 2.0/inv_dVrm

                dxBzm = (U(MAGZ_VAR,i  ,j,k)*xCenter(i) &
                     -  U(MAGZ_VAR,i-1,j,k)*abs(xCenter(i-1)))*inv_dVrm
                Jym = - eta_loc*dxBzm

                eta_loc = 0.5*(res_eta(i,j,k)+res_eta(i+1,j,k))
                inv_dVrp = xCenter(i+1)*abs(xCenter(i+1)) - xCenter(i)*abs(xCenter(i))
                inv_dVrp = 2.0/inv_dVrp

                dxBzp = (U(MAGZ_VAR,i+1,j,k)*xCenter(i+1) &
                     -  U(MAGZ_VAR,i,j,k)*abs(xCenter(i)))*inv_dVrp                 
                Jyp = - eta_loc*dxBzp

                res_source(i,j,k) = -(Jyp-Jym)/dx
             enddo !end of i loop
          enddo !end of j loop
       enddo !end of k loop
    endif
#endif

    iskip = 1
    if (NDIM == 1 .and. rangeSwitch .eq. UPDATE_BOUND) iskip = imax-imin
    do k=kmin,kmax
       do j=jmin,jmax
          if (NDIM >= 2) then
             iskip = 1
             if (rangeSwitch == UPDATE_BOUND .and. j > jmin .and. j < jmax) then
                iskip = imax-imin
                if (NDIM == 3) then
                   iskip = 1
                   if (k > kmin .and. k < kmax) then
                      iskip = imax-imin
                   endif
                endif
             endif
          endif
          do i=imin,imax,iskip
#ifdef BDRY_VAR
             if (U(BDRY_VAR,i,j,k) .LE. 0.0) then
#endif
                !! For non-cartesian geometries
                leftFac = 1.
                rghtFac = 1.

                !! Fluxes at each local cell 
                FL(HY_DENS_FLUX:HY_VOLU_FLUX) = xflux(HY_DENS_FLUX:HY_VOLU_FLUX,i,  j,   k   )
                FR(HY_DENS_FLUX:HY_VOLU_FLUX) = xflux(HY_DENS_FLUX:HY_VOLU_FLUX,i+1,j,   k   )
                GL(HY_DENS_FLUX:HY_VOLU_FLUX) = yflux(HY_DENS_FLUX:HY_VOLU_FLUX,i,  j,   k   )*ky
                GR(HY_DENS_FLUX:HY_VOLU_FLUX) = yflux(HY_DENS_FLUX:HY_VOLU_FLUX,i,  j+ky,k   )*ky
                HL(HY_DENS_FLUX:HY_VOLU_FLUX) = zflux(HY_DENS_FLUX:HY_VOLU_FLUX,i,  j,   k   )*kz
                HR(HY_DENS_FLUX:HY_VOLU_FLUX) = zflux(HY_DENS_FLUX:HY_VOLU_FLUX,i,  j,   k+kz)*kz
#if (NSPECIES+NMASS_SCALARS) > 0
                do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
                   isph= ispu-NPROP_VARS
                   FL(HY_END_FLUX+isph) = xflux(HY_END_FLUX+isph,i,  j,   k   )
                   FR(HY_END_FLUX+isph) = xflux(HY_END_FLUX+isph,i+1,j,   k   )
                   GL(HY_END_FLUX+isph) = yflux(HY_END_FLUX+isph,i,  j,   k   )*ky
                   GR(HY_END_FLUX+isph) = yflux(HY_END_FLUX+isph,i,  j+ky,k   )*ky
                   HL(HY_END_FLUX+isph) = zflux(HY_END_FLUX+isph,i,  j,   k   )*kz
                   HR(HY_END_FLUX+isph) = zflux(HY_END_FLUX+isph,i,  j,   k+kz)*kz
                enddo
#endif

                Sphys = 0.0
                Sphys(HY_XMOM)                  = (FL(HY_P_FLUX) - FR(HY_P_FLUX)) *(1.0-fP)/ dx

                if (hy_geometry /= CARTESIAN) then
                   select case(hy_geometry) ! First, select whether y or z is phi-direction
                   case(CYLINDRICAL)
                      MOM_PHI = HY_ZMOM
                      MOM_PHI_FLUX = HY_ZMOM_FLUX
                      velPhiVar    = VELZ_VAR
                      MOM_ZI       = HY_YMOM
                      MOM_ZI_FLUX  = HY_YMOM_FLUX
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                      MAG_PHI      = HY_MAGZ
                      MAG_PHI_FLUX = HY_MAGZ_FLUX
                      MAG_ZI       = HY_MAGY
                      MAG_ZI_FLUX  = HY_MAGY_FLUX
#endif
                      dz = xCenter(i) * del(DIR_Z)
                      alpha = 1.

                   case(POLAR)
                      MOM_PHI      = HY_YMOM
                      MOM_PHI_FLUX = HY_YMOM_FLUX
                      velPhiVar    = VELY_VAR
                      MOM_ZI       = HY_ZMOM
                      MOM_ZI_FLUX  = HY_ZMOM_FLUX
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                      MAG_PHI      = HY_MAGY
                      MAG_PHI_FLUX = HY_MAGY_FLUX
                      MAG_ZI       = HY_MAGZ
                      MAG_ZI_FLUX  = HY_MAGZ_FLUX
#endif
                      dy = xCenter(i) * del(DIR_Y)
                      alpha = 1.

                   case(SPHERICAL)
                      MOM_PHI      = HY_ZMOM
                      MOM_PHI_FLUX = HY_ZMOM_FLUX
                      velPhiVar    = VELZ_VAR
                      MOM_THT      = HY_YMOM
                      MOM_THT_FLUX = HY_YMOM_FLUX

                      dx_sph = (xRight(i)**3 - xLeft(i)**3) / (3.*xCenter(i)**2)
                      dy     = xCenter(i) * del(DIR_Y)
                      if (NDIM == 3) dz     = xCenter(i) * del(DIR_Z) * sin(yCenter(j))
                      alpha  = 2.
                   end select

                   leftFac = faceAreas(i  ,j,k)*dx/cellVolumes(i,j,k)
                   rghtFac = faceAreas(i+1,j,k)*dx/cellVolumes(i,j,k)


                   if (hy_geometry == CYLINDRICAL) then
#ifndef FLASH_USM_MHD
                      FR = FR*faceAreas(i+1,j,k)
                      FL = FL*faceAreas(i  ,j,k)
                      dxv = cellVolumes(i,j,k)
#endif

#ifdef FLASH_USM_MHD                 
                      FR(HY_DENS_FLUX:HY_MAGY_FLUX) = FR(HY_DENS_FLUX:HY_MAGY_FLUX)*rghtFac
                      FL(HY_DENS_FLUX:HY_MAGY_FLUX) = FL(HY_DENS_FLUX:HY_MAGY_FLUX)*leftFac
                    !!skip MAGZ flux will be treated after this
                      FR(HY_EINT_FLUX:HY_VOLU_FLUX) = FR(HY_EINT_FLUX:HY_VOLU_FLUX)*rghtFac
                      FL(HY_EINT_FLUX:HY_VOLU_FLUX) = FL(HY_EINT_FLUX:HY_VOLU_FLUX)*leftFac
#  if (NSPECIES+NMASS_SCALARS) > 0
                      if (hy_fullSpecMsFluxHandling) then
                         do ispu = SPECIES_BEGIN, MASS_SCALARS_END
                            isph= ispu-NPROP_VARS
                            FR(HY_END_FLUX+isph) = FR(HY_END_FLUX+isph)*rghtFac
                            FL(HY_END_FLUX+isph) = FL(HY_END_FLUX+isph)*leftFac
                         enddo
                      end if
#  endif
#endif

                      !! Angular momentum conservative form
                      if (hy_conserveAngMom) then
                         FR(MOM_PHI_FLUX) = FR(MOM_PHI_FLUX)*rghtFac
                         FL(MOM_PHI_FLUX) = FL(MOM_PHI_FLUX)*leftFac
                      endif


                      !! Angular field conservative form
                      !! This is used when the flag is on.
                      !! This formalism does not require the
                      !! source term in the induction equation
                      !! of Bphi, i.e. Sgeo(MAG_PHI) = 0
#ifdef FLASH_USM_MHD
                      if (.not. hy_conserveAngField) then
                         FR(MAG_PHI_FLUX) = FR(MAG_PHI_FLUX)*rghtFac
                         FL(MAG_PHI_FLUX) = FL(MAG_PHI_FLUX)*leftFac
                      endif
#endif
                   else         ! for SPHERICAL or POLAR geometry
                      FR = FR*faceAreas(i+1,j,k)
                      FL = FL*faceAreas(i  ,j,k)
                      dxv = cellVolumes(i,j,k)
                   endif

                endif !end of non-Cartesian support

                if (NDIM .GE. 2) Sphys(HY_YMOM) = (GL(HY_P_FLUX) - GR(HY_P_FLUX)) *(1.0-fP)/ dy
                if (NDIM .EQ. 3) Sphys(HY_ZMOM) = (HL(HY_P_FLUX) - HR(HY_P_FLUX)) *(1.0-fP)/ dz

#if NDIM > 1
                if (hy_geometry == SPHERICAL) then
                   GR = GR*faceAreasY(i,j+1,k)
                   GL = GL*faceAreasY(i,j  ,k)
                   dyv = cellVolumes(i,j,k)
                else
                   dyv = dy
                end if
#endif

                U0(HY_DENS) = U(DENS_VAR,i,j,k)                                      !density
                if (hy_geometry /= CARTESIAN) then
#ifdef HALFTIME_DENS_FOR_SGEO_TERMS_AS_IN_SPLIT_HYDRO
                   call updateConservedDens&
                     (U0(HY_DENS),densNew,&
                      FL(HY_DENS_FLUX),&
                      FR(HY_DENS_FLUX),&
                      GL(HY_DENS_FLUX),&
                      GR(HY_DENS_FLUX),&
                      HL(HY_DENS_FLUX),&
                      HR(HY_DENS_FLUX),&
                      dxv,dyv,dz,dt)
                   densNph = 0.5 * (U0(HY_DENS) + densNew) ! average of old and new dens
#else
                   densNph = U0(HY_DENS)
#endif
                end if

                !! radially-averaged pressure at n+1/2 (via characteristic tracing, NOT Riemann solver)
                presStar =  scrch_Ptr(HY_VAR2_SCRATCHCTR_VAR,i,j,k)

                Qohm = 0.0
              
#ifdef FLASH_USM_MHD
                if (hy_useMagneticResistivity) then 
                   speciesArr => U(SPECIES_BEGIN:SPECIES_END,i,j,k)
                   call MagneticResistivity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                        speciesArr,res_eta(i,j,k))
                   call hy_uhd_addOhmicHeating(blockID,blkLimits,i,j,k,Qohm,res_eta(i,j,k))
                   Qohm = Qohm*U(DENS_VAR,i,j,k)
                endif
#endif


                if (hy_useAuxEintEqn) then
                   !! Update intenal energy rho*eint
                   IntEner  = U(DENS_VAR,i,j,k)*U(EINT_VAR,i,j,k)
                   !! Note: scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k) holds the volume-averaged pressures at n+1/2
                   tempPres = scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k)

                   call updateInternalEnergy&
                        (IntEner,tempPres,&
                         FL(HY_EINT_FLUX:HY_VOLU_FLUX),&
                         FR(HY_EINT_FLUX:HY_VOLU_FLUX),&
                         GL(HY_EINT_FLUX:HY_VOLU_FLUX),&
                         GR(HY_EINT_FLUX:HY_VOLU_FLUX),&
                         HL(HY_EINT_FLUX:HY_VOLU_FLUX),&
                         HR(HY_EINT_FLUX:HY_VOLU_FLUX),&
                         dxv,dyv,dz,dt,Qohm)
                endif

#ifdef FLASH_UHD_3T
                !! STORE THE OLD DENSITY AND INTERNAL ENERGY FOR 3T UPDATE
                scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,i,j,k) = U(DENS_VAR,i,j,k)
                scrch_Ptr(HY_VAR2_SCRATCHCTR_VAR,i,j,k) = U(EINT_VAR,i,j,k)

#ifdef FLASH_USM_MHD
                scrch_Ptr(HY_XN06_SCRATCHCTR_VAR,i,j,k) = Qohm*dt

                if (hy_hallVelocity) then 
                   !define abar zbar, we store these in scratch for 3T MHD update later on
                   call Eos_getAbarZbar(solnVec=U(:,i,j,k),abar=abar,zbar=zbar) 
                   scrch_Ptr(HY_XN04_SCRATCHCTR_VAR,i,j,k) = abar
                   scrch_Ptr(HY_XN05_SCRATCHCTR_VAR,i,j,k) = zbar

                   ! correct the energy fluxes with current terms
                   !! Note: hy_uhd_getCurrent sets Jp and Jm for many cells, unless called with mode_switch=4.
                   call hy_uhd_getCurrents(blockID, rangeSwitch, blkLimits,datasize, del, Jp, Jm, 4,&
                                           scrch_Ptr,&
                                           i, j, k)
                   Sphys(HY_ENER) = ( Jp(1,i,j,k) - Jm(1,i,j,k) ) / dx
                   if (NDIM > 1) then
                      Sphys(HY_ENER) = Sphys(HY_ENER) + ( Jp(2,i,j,k) - Jm(2,i,j,k) ) / dy
                      if (NDIM == 3) then 
                         HL(HY_ENER_FLUX) = HL(HY_ENER_FLUX) - Jm(3,i,j,k)
                         HR(HY_ENER_FLUX) = HR(HY_ENER_FLUX) - Jp(3,i,j,k)
                      endif
                   endif
                endif
#endif
#endif

                !! Prepare to update conserved quantities
                U0(HY_XMOM:HY_ZMOM) = U(VELX_VAR:VELZ_VAR,i,j,k)*U(DENS_VAR,i,j,k)   !momenta
                U0(HY_ENER) = U(DENS_VAR,i,j,k)*U(ENER_VAR,i,j,k)                    !total gas energy 
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                U0(HY_MAGX:HY_MAGZ) = U(MAGX_VAR:MAGZ_VAR,i,j,k)                     !magnetic fields
                U0(HY_ENER) = U0(HY_ENER)+0.5*dot_product(U(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                          U(MAGX_VAR:MAGZ_VAR,i,j,k))!total plasma energy
#endif
#ifdef FLASH_UGLM_MHD
                U0(HY_GLMP) = U(GLMP_VAR,i,j,k)
#endif

                Sgeo = 0.
                if (hy_geometry /= CARTESIAN) then
                   !! Calculate geometrical source terms.  See S&O 75.
                   !! Advance density and phi-momentum to n+1/2 via finite volume update        
                   densStar = U0(HY_DENS)
                   xvel0    = U(VELX_VAR ,i,j,k)
                   phiVel0  = U(velPhiVar,i,j,k)
                   xmomStar = U0(HY_XMOM)
                   pmomStar = U0(MOM_PHI)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                   xmagStar = U0(HY_MAGX)
                   pmagStar = U0(MAG_PHI)
                   zmagStar = U0(MAG_ZI)
#endif 
                   cs = sqrt(U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k))
                   eta = (abs(U(VELX_VAR,i,j,k)) + cs) * dt/dx
                   eta = (1.-eta) / (cs*dt*abs(alpha/xCenter(i)))

                   Sgeo(HY_XMOM) = (densNph*phiVel0*phiVel0 + fP*alpha*presStar) / xCenter(i)!T phi,phi
                   Sgeo(MOM_PHI) = (densNph*phiVel0*xvel0) / xCenter(i)!T phi,r

                   !! take out the source term for conservative form
                   if (hy_geometry == CYLINDRICAL .and. hy_conserveAngMom) then
                      Sgeo(MOM_PHI) = 0.0
                   endif


#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 ! P* is the total Pressure
                   Sgeo(HY_XMOM) = Sgeo(HY_XMOM) - (pmagStar**2) / xCenter(i)
                   Sgeo(MOM_PHI) = Sgeo(MOM_PHI) - pmagStar*xmagStar / xCenter(i)

                   if (hy_geometry == CYLINDRICAL .and. hy_conserveAngMom) then
                      Sgeo(MOM_PHI) = 0.0
                   endif

                   Sgeo(MAG_PHI) = - (phiVel0 * xmagStar - &
                                     pmagStar * xvel0) / xCenter(i) !O phi,r

                   if (hy_geometry == CYLINDRICAL .and. &
                       hy_conserveAngField) then
                     
                      Sgeo(MAG_PHI) = 0.0 !!- ((pmomStar/densStar) * xmagStar - &
                                   !!pmagStar * (xmomStar/densStar)) / xCenter(i) +&
                                   !!(FR(MAG_PHI_FLUX)-FL(MAG_PHI_FLUX))/dx -&
                                   !!(FR(MAG_PHI_FLUX)*((xCenter(i)+0.5*dx)/xCenter(i))-&
                                   !!FL(MAG_PHI_FLUX)*((xCenter(i)-0.5*dx)/xCenter(i)))/dx
                                   !! !if (rghtFac>0.0 .and. leftFac==0.0) print*,Sgeo(MAG_PHI)
                   endif
                   if (hy_geometry == CYLINDRICAL .and. &
                       hy_useMagneticResistivity) then
                   
                      if (.not. hy_conserveAngField) then
                         Sgeo(MAG_PHI) = Sgeo(MAG_PHI) + res_source(i,j,k) 
                      endif
                   endif
#endif
                   Sgeo(MOM_PHI) = - Sgeo(MOM_PHI)


                   if (hy_geometry == SPHERICAL) then
!!$                      tmomStar = U0(MOM_THT) 
                      thtVel0  = U(VELY_VAR,i,j,k)
                      Sgeo(HY_XMOM) = Sgeo(HY_XMOM) + densNph*thtVel0*thtVel0 / xCenter(i)
                      Sgeo(HY_XMOM) = Sgeo(HY_XMOM)*dx/dx_sph
#if NDIM > 1
                      Sgeo(MOM_THT) = (densNph*phiVel0*phiVel0 + 0*fP*alpha*presStar) / xCenter(i)
                      Sgeo(MOM_THT) = Sgeo(MOM_THT) * cos(yCenter(j))/sin(yCenter(j))
                      Sgeo(MOM_THT) = Sgeo(MOM_THT) - (densNph*thtVel0*xvel0) / xCenter(i)
                      Sgeo(MOM_THT) = Sgeo(MOM_THT)*dx/dx_sph
#endif
                   endif
                endif

                !! Now really update conserved quantities
                call updateConservedVariable&
                     (U0(HY_DENS:HY_DENS+HY_VARINUM-1),&
                      FL(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1),&
                      FR(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1),&
                      GL(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1),&
                      GR(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1),&
                      HL(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1),&
                      HR(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1),&
                      gravX(i,j,k),gravY(i,j,k),gravZ(i,j,k),dxv,dyv,dz,dt,Sgeo,Sphys)

#ifdef DEBUG_HYDRO_POSITIVITY
                if (U0(HY_DENS)<hy_smalldens) then
                   print*,'Low DENS',U(DENS_VAR,i,j,k),'->',U0(HY_DENS),',X=',xCenter(i),',i,j=',i,j,&
                            ' in Block',blockID,'@',hy_meshMe
                end if
#endif


                U(DENS_VAR,i,j,k) = max(U0(HY_DENS),hy_smalldens)                    !density
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                U(MAGX_VAR:MAGZ_VAR,i,j,k) = U0(HY_MAGX:HY_MAGZ)                     !magnetic fields
#endif
#ifdef FLASH_UGLM_MHD
                U(GLMP_VAR,i,j,k) = U0(HY_GLMP)
#endif
#ifdef DEBUG_HYDRO_POSITIVITY
                if (U0(HY_ENER)/U(DENS_VAR,i,j,k)<hy_smallE) then
                   print*,'Low ENER',U(ENER_VAR,i,j,k),'->',U0(HY_ENER)/U(DENS_VAR,i,j,k),',X=',xCenter(i),',i,j=',i,j,&
                            ' in Block',blockID,'@',hy_meshMe
                end if
#endif
                U(ENER_VAR,i,j,k) = U0(HY_ENER)                                      !total plasma energy
                !! We will update velocity fields after species & mass scalar update

#if (NSPECIES+NMASS_SCALARS) > 0
                newDens = U(DENS_VAR,i,j,k)
                if (hy_fullSpecMsFluxHandling) then
                   if (rangeSwitch==UPDATE_ALL_SPECMS_BOUND) then
                      specialForInterior = (i > imin .and. i < imax)
                      if (NDIM > 1) specialForInterior = specialForInterior .AND. (j > jmin .and. j < jmax)
                      if (NDIM > 2) specialForInterior = specialForInterior .AND. (k > kmin .and. k < kmax)
                   else
                      specialForInterior = .FALSE.
                   end if
                   sumSpecies = 0.
                   do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
                      if (specialForInterior) then
                         U(ispu,i,j,k) = U(ispu,i,j,k) / newDens
                      else
                         isph= ispu-NPROP_VARS
                         call updateSpeciesMassScalar&
                                 (U(ispu,i,j,k),Uold(HY_DENS,i,j,k),newDens,&
                                 FL(HY_END_FLUX+isph),FR(HY_END_FLUX+isph),&
                                 GL(HY_END_FLUX+isph),GR(HY_END_FLUX+isph),&
                                 HL(HY_END_FLUX+isph),HR(HY_END_FLUX+isph),&
                                 dxv,dyv,dz,dt)
                      end if

                      !! Conserving mass fractions
                      if (ispu <= SPECIES_END) then
                         sumSpecies = sumSpecies + U(ispu,i,j,k)
                      endif

                   enddo

                !! Conserving mass fractions: They will add up to 1 after this.
                   do ispu = SPECIES_BEGIN, SPECIES_END
                      U(ispu,i,j,k) = U(ispu,i,j,k)/sumSpecies
                   enddo
                else
                !! Note that the velocity fields here are old velocities at time step n, not n+1
                   call hy_uhd_updateSpeciesMassScalar&
                     (hy_order,newDens,&
                      SpOld(1:hy_numXN,i-3:i+3,j-3*ky:j+3*ky,k-3*kz:k+3*kz),&
                      Uold(1:6,        i-3:i+3,j-3*ky:j+3*ky,k-3*kz:k+3*kz),&
                      FL(HY_DENS_FLUX),FR(HY_DENS_FLUX),&
                      GL(HY_DENS_FLUX),GR(HY_DENS_FLUX),&
                      HL(HY_DENS_FLUX),HR(HY_DENS_FLUX),&
                      dxv,dyv,dz,dt, U(SPECIES_BEGIN:MASS_SCALARS_END,i,j,k))
                end if
#endif /*if (NSPECIES+NMASS_SCALARS) > 0*/

                !! Update velocity fields after species & mass scalar update
                !! DL - why should it be AFTER??? Ah... ok, in the old way, 
                !!      it still uses velocity fields for mass scalar and species update
                U(VELX_VAR:VELZ_VAR,i,j,k) = U0(HY_XMOM:HY_ZMOM)/U(DENS_VAR,i,j,k)        !velocities

                U(ENER_VAR,i,j,k) = U0(HY_ENER)/U(DENS_VAR,i,j,k) !total plasma energy, now in mass-specific form

                if (hy_useAuxEintEqn) then
                   IntEner = IntEner / U(DENS_VAR,i,j,k)
#ifdef EINT_VAR
#ifdef DEBUG_HYDRO_POSITIVITY
                   if (IntEner<hy_smallE) then
                      print*,'Low EINT',U(EINT_VAR,i,j,k),'->',IntEner,',X=',xCenter(i),',i,j=',i,j,&
                           ' in Block',blockID,'@',hy_meshMe
                   end if
#endif
                   U(EINT_VAR,i,j,k) = IntEner
#endif
                end if

#ifdef FLASH_UHD_3T
#ifdef FLASH_UHD_HYDRO
                !! Perform energy updates here for 3T.  For 1T energy is updated in hy_uhd_energyFix.
                !! We only perform the energy updates here for pure hydro only.
                !! MHD still needs to update magnetic fields after the current cell-centered
                !! variable updates. Therefore, the energy updates should be done after
                !! updating magnetic fields.
                ekin = .5*dot_product(U(VELX_VAR:VELZ_VAR,i,j,k),U(VELX_VAR:VELZ_VAR,i,j,k))

                eint = U(ENER_VAR,i,j,k)-ekin

                if (.not. hy_useAuxEintEqn .or. eint > hy_eswitch*ekin) then
                   newEint = max(hy_smallE,eint)
                else
                   newEint = max(hy_smallE,IntEner)
                endif
                !! Store specific gas energy ener = ekin + eint
                U(ENER_VAR,i,j,k) = newEint + ekin
#ifdef EINT_VAR
                U(EINT_VAR,i,j,k) = newEint
#endif
#endif /*ifdef FLASH_UHD_HYDRO*/
#endif /*ifdef FLASH_UHD_3T*/

#ifdef BDRY_VAR
             endif
#endif

          enddo !end of i loop
       enddo !end of j loop
    enddo !end of k loop

    !! ---------------------------------------------------------------
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
    if (hy_forceHydroLimit) then
       U(MAGX_VAR:MAGZ_VAR,:,:,:) = 0.
    endif
#endif

#if (NSPECIES+NMASS_SCALARS) > 0
  ! Renormalize or limit abundances
    if (hy_irenorm == 1) then
       call Grid_renormAbundance(blockID,blkLimits,U)
    else
       call Grid_limitAbundance(blkLimits,U)
    endif
#endif  

    ! Release block pointers
    call Grid_releaseBlkPtr(blockID,U,CENTER)



  CONTAINS

    !  The rangeSwitch==UPDATE_SPECMS_INTERIOR case is farmed out to this simpler routine.
    subroutine unsplitUpdateSpecMs

#if (NSPECIES+NMASS_SCALARS) > 0


      dx = del(DIR_X)
      dy = 1.
      dz = 1.
      if (NDIM >= 2) then
         dy = del(DIR_Y)
         if (NDIM == 3) then
            dz = del(DIR_Z)
         endif
      endif
      dxv = dx; dyv = dy

    !! Set region to update for update mode UPDATE_SPECMS_INTERIOR
      imin  = blkLimits(LOW, IAXIS)+1
      imax  = blkLimits(HIGH,IAXIS)-1
      if (NDIM >= 2) then
         jmin  = blkLimits(LOW, JAXIS)+1
         jmax  = blkLimits(HIGH,JAXIS)-1
         if (NDIM == 3) then
            kmin  = blkLimits(LOW, KAXIS)+1
            kmax  = blkLimits(HIGH,KAXIS)-1
         else
            kmin  = 1
            kmax  = 1
         endif
      else
         jmin  = 1
         jmax  = 1
      endif


    ! Get block pointers
    call Grid_getBlkPtr(blockID,U,CENTER)

    if (hy_geometry /= CARTESIAN) then
       faceAreas = 0.
       call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, &
            (/blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS)/), &
            faceAreas(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),  &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)), &
            (/isize+1, jsize, ksize/) )
#if NDIM > 1
       if (hy_geometry == SPHERICAL) then
          call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
            (/blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS)/), &
            faceAreasY(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,  &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)), &
            (/isize, jsize+1, ksize/) )
       end if
#endif

       cellVolumes = 0.
       call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR, &
            (/blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS)/), &
            cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)), &
            (/isize, jsize, ksize/) )
    endif

    

    if (.NOT. hy_fullSpecMsFluxHandling) then
       do ispu =  SPECIES_BEGIN, MASS_SCALARS_END
          isph= ispu-NPROP_VARS
          SpOld(isph,:,:,:) = U(ispu,:,:,:)
       enddo
       Uold(1,  :,:,:) = U(DENS_VAR,:,:,:)
       Uold(2,  :,:,:) = U(PRES_VAR,:,:,:)
       Uold(3:5,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)
       Uold(6,  :,:,:) = U(GAME_VAR,:,:,:)
    end if

    if (hy_geometry /= CARTESIAN) then
       call Grid_getCellCoords(IAXIS,blockID, CENTER,    .true.,xCenter, dataSize(IAXIS))
       call Grid_getCellCoords(IAXIS,blockID, LEFT_EDGE, .true.,xLeft,   dataSize(IAXIS))
       call Grid_getCellCoords(IAXIS,blockID, RIGHT_EDGE,.true.,xRight,  dataSize(IAXIS))
       if (NDIM == 3 .AND. hy_geometry == SPHERICAL) then
          call Grid_getCellCoords(JAXIS,blockID, CENTER,.false.,yCenter, size(yCenter))
       end if
    endif


    ! define dimension dependent switches
    ky=0
    kz=0

    if (NDIM > 1) then
       ky=1
       if (NDIM > 2) then
          kz=1
       endif
    endif

    do k=kmin,kmax
       do j=jmin,jmax
          do i=imin,imax
#ifdef BDRY_VAR
             if (U(BDRY_VAR,i,j,k) .LE. 0.0) then
#endif
                !! For non-cartesian geometries
                leftFac = 1.
                rghtFac = 1.

                !! Fluxes at each local cell 
                do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
                   isph= ispu-NPROP_VARS
                   FL(HY_END_FLUX+isph) = xflux(HY_END_FLUX+isph,i,  j,   k   )
                   FR(HY_END_FLUX+isph) = xflux(HY_END_FLUX+isph,i+1,j,   k   )
                   GL(HY_END_FLUX+isph) = yflux(HY_END_FLUX+isph,i,  j,   k   )*ky
                   GR(HY_END_FLUX+isph) = yflux(HY_END_FLUX+isph,i,  j+ky,k   )*ky
                   HL(HY_END_FLUX+isph) = zflux(HY_END_FLUX+isph,i,  j,   k   )*kz
                   HR(HY_END_FLUX+isph) = zflux(HY_END_FLUX+isph,i,  j,   k+kz)*kz
                enddo



                if (hy_geometry /= CARTESIAN) then
                   if (NDIM > 1) then
                      select case(hy_geometry) ! First, select whether y or z is phi-direction
                      case(CYLINDRICAL)
                         if (NDIM == 3) dz = xCenter(i) * del(DIR_Z)

                      case(POLAR)
                         dy = xCenter(i) * del(DIR_Y)

                      case(SPHERICAL)
                         dy     = xCenter(i) * del(DIR_Y)
                         if (NDIM == 3) dz     = xCenter(i) * del(DIR_Z) * sin(yCenter(j))
                      end select
                   end if

                   leftFac = faceAreas(i  ,j,k)
                   rghtFac = faceAreas(i+1,j,k)
                   dxv = cellVolumes(i,j,k)

                   if (hy_fullSpecMsFluxHandling) then
                      do ispu = SPECIES_BEGIN, MASS_SCALARS_END
                         isph= ispu-NPROP_VARS
                         FR(HY_END_FLUX+isph) = FR(HY_END_FLUX+isph)*rghtFac
                         FL(HY_END_FLUX+isph) = FL(HY_END_FLUX+isph)*leftFac
                      enddo
                   else
                      FR(HY_DENS_FLUX) = FR(HY_DENS_FLUX)*rghtFac
                      FL(HY_DENS_FLUX) = FL(HY_DENS_FLUX)*leftFac
                   end if

                   dyv = dy
#if NDIM > 1
                   if (hy_geometry == SPHERICAL) then
                      leftFac = faceAreasY(i,j  ,k)
                      rghtFac = faceAreasY(i,j+1,k)
                      dyv = cellVolumes(i,j,k)

                      if (hy_fullSpecMsFluxHandling) then
                         do ispu = SPECIES_BEGIN, MASS_SCALARS_END
                            isph= ispu-NPROP_VARS
                            GR(HY_END_FLUX+isph) = GR(HY_END_FLUX+isph)*rghtFac
                            GL(HY_END_FLUX+isph) = GL(HY_END_FLUX+isph)*leftFac
                         enddo
                      else
                         GR(HY_DENS_FLUX) = GR(HY_DENS_FLUX)*rghtFac
                         GL(HY_DENS_FLUX) = GL(HY_DENS_FLUX)*leftFac
                      end if
                   end if
#endif

                endif !end of non-Cartesian support


                if (hy_fullSpecMsFluxHandling) then
                   do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
                      isph= ispu-NPROP_VARS
                      call updateSpeciesMassScalar&
                                 (U(ispu,i,j,k),U(DENS_VAR,i,j,k), 1.0,&
                                 FL(HY_END_FLUX+isph),FR(HY_END_FLUX+isph),&
                                 GL(HY_END_FLUX+isph),GR(HY_END_FLUX+isph),&
                                 HL(HY_END_FLUX+isph),HR(HY_END_FLUX+isph),&
                                 dxv,dyv,dz,dt)
                   enddo

                else
                !! Note that the velocity fields here are old velocities at time step n, not n+1
                   call hy_uhd_updateSpeciesMassScalar&
                     (hy_order,0.0,&
                      SpOld(1:hy_numXN,i-3:i+3,j-3*ky:j+3*ky,k-3*kz:k+3*kz),&
                      Uold(1:6,        i-3:i+3,j-3*ky:j+3*ky,k-3*kz:k+3*kz),&
                      FL(HY_DENS_FLUX),FR(HY_DENS_FLUX),&
                      GL(HY_DENS_FLUX),GR(HY_DENS_FLUX),&
                      HL(HY_DENS_FLUX),HR(HY_DENS_FLUX),&
                      dxv,dyv,dz,dt, U(SPECIES_BEGIN:MASS_SCALARS_END,i,j,k))
                end if


#ifdef BDRY_VAR
             endif
#endif

          enddo !end of i loop
       enddo !end of j loop
    enddo !end of k loop

    !! ---------------------------------------------------------------
    ! Release block pointers
    call Grid_releaseBlkPtr(blockID,U,CENTER)
#endif

  End Subroutine unsplitUpdateSpecMs







  End Subroutine hy_uhd_unsplitUpdate



!! ==================================================================
Subroutine updateConservedVariable(Ul,FL,FR,GL,GR,HL,HR,  &
                                      gravX, gravY, gravZ,&
                                      dx,dy,dz,dt,Sgeo,Sphys)
  use Hydro_data, ONLY : hy_useGravity
  implicit none
  real, dimension(HY_VARINUM), intent(INOUT) :: Ul
  real, dimension(HY_VARINUM), intent(IN) :: Sgeo,Sphys
  real, dimension(HY_VARINUM), intent(IN) :: FL,FR,GL,GR,HL,HR
  real, intent(IN) :: gravX,gravY,gravZ
  real, intent(IN) :: dx,dy,dz,dt
  real, dimension(3) :: momentaOld
  real :: densOld

  !! Store old states at n
  densOld = Ul(HY_DENS)
  momentaOld(1:3) = Ul(HY_XMOM:HY_ZMOM)

  !! Update conservative variables from n to n+1 step
  Ul(HY_DENS:HY_DENS+HY_VARINUM-1)=Ul(HY_DENS:HY_DENS+HY_VARINUM-1)&
       -dt/dx*( FR(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1)&
               -FL(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1))
  
  if (NDIM > 1) then
  Ul(HY_DENS:HY_DENS+HY_VARINUM-1)=Ul(HY_DENS:HY_DENS+HY_VARINUM-1)&
       -dt/dy*( GR(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1)&
               -GL(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1))

  if (NDIM > 2) then
  Ul(HY_DENS:HY_DENS+HY_VARINUM-1)=Ul(HY_DENS:HY_DENS+HY_VARINUM-1)& 
       -dt/dz*( HR(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1)&
               -HL(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1))
  endif

  endif


  !! Include geometric and physical source terms
  Ul = Ul + dt*(Sgeo+Sphys)

  if (hy_useGravity) then
     !! Note: 1. Extrapolated gravity case
     !!          gravX,Y,Z - gravity at n
     !!------------------------------------------

     !! Update total energy (rho*E) with gravity source
     !! The following is new gravity formulation.  Still needs to be corrected by n+1 gravity!!!
     Ul(HY_XMOM:HY_ZMOM) = Ul(HY_XMOM:HY_ZMOM)&
          + 0.5*dt*densOld*(/gravX,gravY,gravZ/)

     Ul(HY_ENER) = Ul(HY_ENER) &
          + 0.5*dt*dot_product(momentaOld(1:3),(/gravX,gravY,gravZ/))

  endif


End Subroutine updateConservedVariable

!! ==================================================================
Subroutine updateConservedDens(densOld,densOut,FL,FR,GL,GR,HL,HR,  &
                                      dx,dy,dz,dt)
  use Hydro_data, ONLY : hy_useGravity
  implicit none
  real, intent(IN)  :: densOld
  real, intent(OUT) :: densOut
  real, intent(IN) :: FL,FR,GL,GR,HL,HR
  real, intent(IN) :: dx,dy,dz,dt

  real :: densNew

  !! Copy old density at n
  densNew = densOld

  !! Update conservative variables from n to n+1 step
  densNew=densNew&
       -dt/dx*( FR-FL)
  
  if (NDIM > 1) then
     densNew=densNew&
       -dt/dy*( GR-GL)

     if (NDIM > 2) then
        densNew=densNew& 
             -dt/dz*( HR-HL)
     endif

  endif

  densOut = densNew
End Subroutine updateConservedDens


!! ==================================================================
Subroutine updateInternalEnergy(eint,pres,FL,FR,GL,GR,HL,HR,dx,dy,dz,dt,Qohm)

  implicit none

  real, intent(INOUT) :: eint ! (=rho*eint)
  real, dimension(2), intent(IN) :: FL,FR,GL,GR,HL,HR
  real, intent(IN) :: pres,dx,dy,dz,dt,Qohm

  eint = eint + dt/dx*(FL(1)-FR(1) + pres*(FL(2)-FR(2)))

  if (NDIM > 1) then
  eint = eint + dt/dy*(GL(1)-GR(1) + pres*(GL(2)-GR(2)))

  if (NDIM > 2) then
  eint = eint + dt/dz*(HL(1)-HR(1) + pres*(HL(2)-HR(2)))
  endif
  endif
  eint = eint + dt*Qohm 

End Subroutine updateInternalEnergy
!! ==================================================================
Subroutine updateSpeciesMassScalar(Spc,densOld,densNew,FL,FR,GL,GR,HL,HR,dx,dy,dz,dt)

  implicit none

  real, intent(INOUT) :: Spc
  real, intent(IN)    :: densOld,densNew
  real, intent(IN)    :: FL,FR,GL,GR,HL,HR
  real, intent(IN)    :: dx,dy,dz,dt

#if NDIM == 1
  Spc = (Spc*densOld + dt*(FL-FR)/dx)/densNew
#elif NDIM == 2
  Spc = (Spc*densOld + dt*((FL-FR)/dx+(GL-GR)/dy))/densNew
#else
  Spc = (Spc*densOld + dt*((FL-FR)/dx+(GL-GR)/dy+(HL-HR)/dz))/densNew
#endif


End Subroutine updateSpeciesMassScalar
!! ==================================================================
