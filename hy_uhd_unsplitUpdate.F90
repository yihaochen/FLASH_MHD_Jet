!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_unsplitUpdate
!!
!! NAME
!!
!!  hy_uhd_unsplitUpdate
!!
!! SYNOPSIS
!!
!!  hy_uhd_unsplitUpdate( integer(IN) :: blockID,
!!                        integer(IN) :: range_switch,
!!                        real(IN)    :: dt,
!!                        real(IN)    :: del(MDM),
!!                        integer(IN) :: dataSize(3),
!!                        integer(IN) :: blkLimits(2,MDIM),
!!                        integer(IN) :: blkLimitsGC(2,MDIM),
!!                        real(IN)    :: xflux(:,:,:,:), 
!!                        real(IN)    :: yflux(:,:,:,:), 
!!                        real(IN)    :: zflux(:,:,:,:),
!!                        real(IN)    :: gravX(:,:,:),
!!                        real(IN)    :: gravY(:,:,:),
!!                        real(IN)    :: gravZ(:,:,:))
!!
!! ARGUMENTS
!!
!!   blockID      - current block ID
!!   range_switch - switch for selective updates on AMR (only valid for hydro and not for MHD) 
!!   dt           - timestep
!!   del          - deltas in {x,y,z} directions
!!   dataSize     - size of the current block
!!   blkLimits    - an array that holds the lower and upper indices of the section
!!                  of block without the guard cells
!!   blkLimitsGC    - an array that holds the lower and upper indices of the section
!!                  of block with the guard cells
!!   xflux,yflux,zflux - cell face centered fluxes at each {=x,y,z} direction
!!   gravX,gravY,gravZ - gravity components in x,y,z directions
!!
!! DESCRIPTION
!!
!!   This routine updates the cell-centered conservative variables and intermediate
!!   internal energy to the next time step using an directionally unsplit scheme.
!!
!!***

!!REORDER(4): U, Uold, SpOld, scrch_Ptr, [xyz]flux

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "UHD.h"

  Subroutine hy_uhd_unsplitUpdate(blockID,range_switch,dt,del,dataSize,blkLimits,&
                                  blkLimitsGC,xflux,yflux,zflux,gravX,gravY,gravZ)

    use Hydro_data,           ONLY : hy_smalldens,hy_order,hy_irenorm,hy_numXN, &
                                     hy_unsplitEosMode, hy_geometry, hy_gcMaskSize, &
                                     hy_gcMask, hy_eswitch, hy_useAuxEintEqn, &
                                     hy_conserveAngMom, hy_smallE
#ifdef FLASH_USM_MHD
  use Hydro_data,             ONLY : hy_mref, hy_E_upwind, hy_hallVelocity, &
                                     hy_useMagneticResistivity, hy_conserveAngField
  use MagneticResistivity_interface, &
                              ONLY : MagneticResistivity
  use hy_uhd_interface,       ONLY : hy_uhd_addOhmicHeating 
#endif
#if defined(FLASH_USM_MHD)  
    use Hydro_data,           ONLY : hy_forceHydroLimit
#endif
    use hy_uhd_interface,     ONLY : hy_uhd_updateSpeciesMassScalar
#ifdef FLASH_UHD_3T
#ifdef FLASH_USM_MHD
  use hy_uhd_interface,       ONLY : hy_uhd_getCurrents
#endif
#endif
    use Grid_interface,       ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
                                     Grid_getCellCoords, Grid_getBlkData
    use Hydro_interface,      ONLY : Hydro_recalibrateEints
    use hy_uhd_slopeLimiters, ONLY : checkMedian, mc, minmod, vanLeer
#ifdef FLASH_UHD_3T
#ifdef FLASH_USM_MHD
    use Eos_interface,        ONLY : Eos_getAbarZbar
#endif
#endif

    implicit none

    !! ---- Arguments ---------------------------------
    integer,intent(IN) :: blockID, range_switch
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
    real, dimension(HY_VARINUM) :: U0
    real, dimension(NFLUXES)    :: FL,FR,GL,GR,HL,HR
    real, pointer, dimension(:,:,:,:) :: U, scrch_Ptr
    real    :: IntEner,tempPres
    integer :: sp, spn

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
!<- ychen 10-2014
    real, dimension(dataSize(IAXIS)) :: xcent
    real, dimension(dataSize(JAXIS)) :: ycent
    real, dimension(dataSize(KAXIS)) :: zcent
! ychen ->

    logical :: gcMask(hy_gcMaskSize)
    real, dimension(HY_VARINUM) :: Sgeo
    real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)) :: faceAreas, cellVolumes
    real, dimension(3, dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)) :: Jp, Jm
    integer :: kx,ky,kz,iskip
    real    :: temp, hdt
    real    :: presStar, densStar, pmomStar, tmomStar, xmomStar
    real    :: pmagStar, xmagStar, zmagStar
    integer :: VEL_PHI, MOM_PHI, MOM_PHI_FLUX, MAG_PHI,  MAG_PHI_FLUX
    integer :: VEL_ZI, MOM_ZI, MOM_ZI_FLUX, MAG_ZI,  MAG_ZI_FLUX
    integer :: VEL_THT, MOM_THT, MOM_THT_FLUX
    real    :: leftFac, rghtFac, dPdr, rvol, alpha, cs, eta
    real    :: ekin, eint
    real    :: abar, zbar
    real    :: Qohm
    integer :: isph,ispu
    real    :: sumSpecies

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

#ifndef FLASH_USM_MHD
  !This is just here so that we can maintain a single OpenMP
  !parallel region for both hydro and mhd applications.
  logical :: hy_E_upwind
#endif

    hdt = 0.5*dt

    iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
    jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
    kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

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

    !! Set regions to update depending on update mode
    iskip = 1
    if (range_switch==UPDATE_INTERIOR) then
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
    call Grid_getBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)

    if (hy_geometry /= CARTESIAN) then
       faceAreas = 0.
       call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, &
            (/blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS)/), &
            faceAreas(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),  &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)), &
            (/isize+1, jsize, ksize/) )

       cellVolumes = 0.
       call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR, &
            (/blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS)/), &
            cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
            blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
            blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)), &
            (/isize, jsize, ksize/) )
    endif

    
#if defined(FLASH_USM_MHD) 
    if (hy_forceHydroLimit) then
       U(MAGX_VAR:MAGZ_VAR,:,:,:) = 0.
    endif
#endif


#if (NSPECIES+NMASS_SCALARS) > 0
    do sp =  SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
       spn= sp-NPROP_VARS
       SpOld(spn,:,:,:) = U(sp,:,:,:)
    enddo
    Uold(1,  :,:,:) = U(DENS_VAR,:,:,:)
    Uold(2,  :,:,:) = U(PRES_VAR,:,:,:)
    Uold(3:5,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)
    Uold(6,  :,:,:) = U(GAME_VAR,:,:,:)
#endif

    if (hy_geometry /= CARTESIAN) then
       call Grid_getCellCoords(IAXIS,blockID, CENTER,    .true.,xCenter, dataSize(IAXIS))
       call Grid_getCellCoords(IAXIS,blockID, LEFT_EDGE, .true.,xLeft,   dataSize(IAXIS))
       call Grid_getCellCoords(IAXIS,blockID, RIGHT_EDGE,.true.,xRight,  dataSize(IAXIS))
    endif

! <- ychen 10-2014
    call Grid_getCellCoords(IAXIS,blockID, CENTER,    .true.,xcent, dataSize(IAXIS))
    call Grid_getCellCoords(JAXIS,blockID, CENTER,    .true.,ycent, dataSize(JAXIS))
    call Grid_getCellCoords(KAXIS,blockID, CENTER,    .true.,zcent, dataSize(KAXIS))
! ychen ->

#ifdef FLASH_USM_MHD
#ifdef FLASH_UHD_3T
    !! STORE THE OLD FIELD FOR CURRENT CALCULATION FOR 3T UPDATE
    scrch_Ptr(XN01_SCRATCH_CENTER_VAR,:,:,:) = U(MAGX_VAR,:,:,:)
    scrch_Ptr(XN02_SCRATCH_CENTER_VAR,:,:,:) = U(MAGY_VAR,:,:,:)
    scrch_Ptr(XN03_SCRATCH_CENTER_VAR,:,:,:) = U(MAGZ_VAR,:,:,:)
#endif
#endif

    ! define dimension dependent switches
    kx=0
    ky=0
    kz=0

    if (NDIM > 1) then
       ky=1
       if (NDIM > 2) then
          kz=1
       endif
    endif

#ifdef FLASH_USM_MHD
    if (hy_E_upwind) kx = 1
#endif


    iskip = 1
    if (NDIM == 1 .and. range_switch .eq. UPDATE_BOUND) iskip = imax-imin
!!! Loop to get magnetic resistivity source correction 
#ifdef FLASH_USM_MHD
    if (hy_geometry == CYLINDRICAL .and. hy_useMagneticResistivity) then
       do k=kmin-kx*kz,kmax+kx*kz
          do j=jmin-kx*ky-2,jmax+kx*ky+2
             if (NDIM >= 2) then
                iskip = 1
                if (range_switch == UPDATE_BOUND .and. j > jmin .and. j < jmax) then
                   iskip = imax-imin
                   if (NDIM == 3) then
                      iskip = 1
                      if (k > kmin .and. k < kmax) then
                         iskip = imax-imin
                      endif
                   endif
                endif
             endif

             do i=imin-kx-2,imax+kx+2,iskip
                !! Get magnetic Resistivity
                speciesArr => U(SPECIES_BEGIN:SPECIES_END,i-1,j,k)
                call MagneticResistivity(U(TEMP_VAR,i-1,j,k),U(DENS_VAR,i-1,j,k),&
                     speciesArr,res_eta(i-1,j,k))
                speciesArr => U(SPECIES_BEGIN:SPECIES_END,i,j,k)
                call MagneticResistivity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                     speciesArr,res_eta(i,j,k))
                speciesArr => U(SPECIES_BEGIN:SPECIES_END,i+1,j,k)
                call MagneticResistivity(U(TEMP_VAR,i+1,j,k),U(DENS_VAR,i+1,j,k),&
                     speciesArr,res_eta(i+1,j,k))                  
                !normalize if needed
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
    if (NDIM == 1 .and. range_switch .eq. UPDATE_BOUND) iskip = imax-imin
    do k=kmin-kx*kz,kmax+kx*kz
       do j=jmin-kx*ky,jmax+kx*ky
          if (NDIM >= 2) then
             iskip = 1
             if (range_switch == UPDATE_BOUND .and. j > jmin .and. j < jmax) then
                iskip = imax-imin
                if (NDIM == 3) then
                   iskip = 1
                   if (k > kmin .and. k < kmax) then
                      iskip = imax-imin
                   endif
                endif
             endif
          endif
          do i=imin-kx,imax+kx,iskip
#ifdef BDRY_VAR
             if (U(BDRY_VAR,i,j,k) < 0.0) then
#endif
                !! For non-cartesian geometries
                leftFac = 1.
                rghtFac = 1.

                !! Fluxes at each local cell 
                FL(HY_DENS_FLUX:HY_PRES_FLUX) = xflux(HY_DENS_FLUX:HY_PRES_FLUX,i,  j,   k   )
                FR(HY_DENS_FLUX:HY_PRES_FLUX) = xflux(HY_DENS_FLUX:HY_PRES_FLUX,i+1,j,   k   )
                GL(HY_DENS_FLUX:HY_PRES_FLUX) = yflux(HY_DENS_FLUX:HY_PRES_FLUX,i,  j,   k   )*ky
                GR(HY_DENS_FLUX:HY_PRES_FLUX) = yflux(HY_DENS_FLUX:HY_PRES_FLUX,i,  j+ky,k   )*ky
                HL(HY_DENS_FLUX:HY_PRES_FLUX) = zflux(HY_DENS_FLUX:HY_PRES_FLUX,i,  j,   k   )*kz
                HR(HY_DENS_FLUX:HY_PRES_FLUX) = zflux(HY_DENS_FLUX:HY_PRES_FLUX,i,  j,   k+kz)*kz
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



                if (hy_geometry /= CARTESIAN) then
                   select case(hy_geometry) ! First, select whether y or z is phi-direction
                   case(CYLINDRICAL)
                      MOM_PHI = HY_ZMOM
                      MOM_PHI_FLUX = HY_ZMOM_FLUX
                      MOM_ZI       = HY_YMOM
                      MOM_ZI_FLUX  = HY_YMOM_FLUX
#if defined(FLASH_USM_MHD) 
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
                      MOM_ZI       = HY_ZMOM
                      MOM_ZI_FLUX  = HY_ZMOM_FLUX
#if defined(FLASH_USM_MHD) 
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
                      MOM_THT      = HY_YMOM
                      MOM_THT_FLUX = HY_YMOM_FLUX

                      dx_sph = (xRight(i)**3 - xLeft(i)**3) / (3.*xCenter(i)**2)
                      dy     = xCenter(i) * del(DIR_Y)
                      dz     = xCenter(i) * del(DIR_Z)
                      alpha  = 2.
                   end select

                   leftFac = faceAreas(i  ,j,k)*dx/cellVolumes(i,j,k)
                   rghtFac = faceAreas(i+1,j,k)*dx/cellVolumes(i,j,k)


                   if (hy_geometry == CYLINDRICAL) then
#ifndef FLASH_USM_MHD
                    FR = FR*rghtFac
                    FL = FL*leftFac
#endif

#ifdef FLASH_USM_MHD                 
                    FR(HY_DENS_FLUX:HY_MAGY_FLUX) = FR(HY_DENS_FLUX:HY_MAGY_FLUX)*rghtFac
                    FL(HY_DENS_FLUX:HY_MAGY_FLUX) = FL(HY_DENS_FLUX:HY_MAGY_FLUX)*leftFac
                    !!skip MAGZ flux will be treated after this
                    FR(HY_EINT_FLUX:HY_PRES_FLUX) = FR(HY_EINT_FLUX:HY_PRES_FLUX)*rghtFac
                    FL(HY_EINT_FLUX:HY_PRES_FLUX) = FL(HY_EINT_FLUX:HY_PRES_FLUX)*leftFac
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
                   endif

                   if (hy_geometry == SPHERICAL) then
                      FR = FR*dx_sph/dx
                      FL = FL*dx_sph/dx
                   endif

                endif !end of non-Cartesian support

                !! radially-averaged pressure at n+1/2 (via charecteristic tracing, NOT Riemann solver)
                presStar =  scrch_Ptr(VAR2_SCRATCH_CENTER_VAR,i,j,k)

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
                   !! Note: scrch_Ptr(VAR1_SCRATCH_CENTER_VAR,i,j,k) holds the volume-averaged pressures at n+1/2
                   tempPres = scrch_Ptr(VAR1_SCRATCH_CENTER_VAR,i,j,k)

                   call updateInternalEnergy&
                        (IntEner,tempPres,&
                         FL(HY_EINT_FLUX:HY_PRES_FLUX),&
                         FR(HY_EINT_FLUX:HY_PRES_FLUX),&
                         GL(HY_EINT_FLUX:HY_PRES_FLUX),&
                         GR(HY_EINT_FLUX:HY_PRES_FLUX),&
                         HL(HY_EINT_FLUX:HY_PRES_FLUX),&
                         HR(HY_EINT_FLUX:HY_PRES_FLUX),&
                         dx,dy,dz,dt,Qohm)
                endif

#ifdef FLASH_UHD_3T
                !! STORE THE OLD DENSITY AND INTERNAL ENERGY FOR 3T UPDATE
                scrch_Ptr(VAR1_SCRATCH_CENTER_VAR,i,j,k) = U(DENS_VAR,i,j,k)
                scrch_Ptr(VAR2_SCRATCH_CENTER_VAR,i,j,k) = U(EINT_VAR,i,j,k)

#ifdef FLASH_USM_MHD
                scrch_Ptr(XN06_SCRATCH_CENTER_VAR,i,j,k) = Qohm*dt

                if (hy_hallVelocity) then 
                   !define abar zbar, we store these in scratch for 3T MHD update later on
                   call Eos_getAbarZbar(solnVec=U(:,i,j,k),abar=abar,zbar=zbar) 
                   scrch_Ptr(XN04_SCRATCH_CENTER_VAR,i,j,k) = abar
                   scrch_Ptr(XN05_SCRATCH_CENTER_VAR,i,j,k) = zbar

                   ! correct the energy fluxes with current terms
                   !! PT: hy_uhd_getCurrent sets Jp and Jm for many cells, should therefore not be called inside an i,j,k loop!
                   call hy_uhd_getCurrents(blockID, range_switch, blkLimits,datasize, del, Jp, Jm, 4, i, j, k)
                   FL(HY_ENER_FLUX) = FL(HY_ENER_FLUX) - Jm(1,i,j,k)
                   FR(HY_ENER_FLUX) = FR(HY_ENER_FLUX) - Jp(1,i,j,k)
                   if (NDIM > 1) then
                      GL(HY_ENER_FLUX) = GL(HY_ENER_FLUX) - Jm(2,i,j,k)
                      GR(HY_ENER_FLUX) = GR(HY_ENER_FLUX) - Jp(2,i,j,k)
                      if (NDIM == 3) then 
                         HL(HY_DENS_FLUX) = HL(HY_DENS_FLUX) - Jm(3,i,j,k)
                         HR(HY_DENS_FLUX) = HR(HY_DENS_FLUX) - Jp(3,i,j,k)
                      endif
                   endif
                endif
#endif
#endif
                if (hy_useAuxEintEqn) U(EINT_VAR,i,j,k) = IntEner ! rho*eint

                !! Update conserved quantities
                U0(HY_DENS) = U(DENS_VAR,i,j,k)                                      !density
                U0(HY_XMOM:HY_ZMOM) = U(VELX_VAR:VELZ_VAR,i,j,k)*U(DENS_VAR,i,j,k)   !momenta
                U0(HY_ENER) = U(DENS_VAR,i,j,k)*U(ENER_VAR,i,j,k)                    !total gas energy 
#if defined(FLASH_USM_MHD) 
                U0(HY_MAGX:HY_MAGZ) = U(MAGX_VAR:MAGZ_VAR,i,j,k)                     !magnetic fields
                U0(HY_ENER) = U0(HY_ENER)+0.5*dot_product(U(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                          U(MAGX_VAR:MAGZ_VAR,i,j,k))!total plasma energy
#endif

                Sgeo = 0.
                if (hy_geometry /= CARTESIAN) then
                   !! Calculate geometrical source terms.  See S&O 75.
                   !! Advance density and phi-momentum to n+1/2 via finite volume update        
                   densStar = U0(HY_DENS)
                   xmomStar = U0(HY_XMOM)
                   pmomStar = U0(MOM_PHI)

#if defined(FLASH_USM_MHD) 
                   xmagStar = U0(HY_MAGX)
                   pmagStar = U0(MAG_PHI)
                   zmagStar = U0(MAG_ZI)
#endif 
                   cs = sqrt(U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k))
                   eta = (abs(U(VELX_VAR,i,j,k)) + cs) * dt/dx
                   eta = (1.-eta) / (cs*dt*abs(alpha/xCenter(i)))

                   Sgeo(HY_XMOM) = (pmomStar**2/densStar + alpha*presStar) / xCenter(i)!T phi,phi
                   Sgeo(MOM_PHI) = (pmomStar*xmomStar/densStar) / xCenter(i)!T phi,r

                   !! take out the source term for conservative form
                   if (hy_geometry == CYLINDRICAL .and. hy_conserveAngMom) then
                      Sgeo(MOM_PHI) = 0.0
                   endif


#if defined(FLASH_USM_MHD) 
                 ! P* is the total Pressure
                 Sgeo(HY_XMOM) = Sgeo(HY_XMOM) - (pmagStar**2) / xCenter(i)
                 Sgeo(MOM_PHI) = Sgeo(MOM_PHI) - pmagStar*xmagStar / xCenter(i)

                 if (hy_geometry == CYLINDRICAL .and. hy_conserveAngMom) then
                   Sgeo(MOM_PHI) = 0.0
                 endif

                   Sgeo(MAG_PHI) = - ((pmomStar/densStar) * xmagStar - &
                                     pmagStar * (xmomStar/densStar)) / xCenter(i) !O phi,r

                 if (hy_geometry == CYLINDRICAL .and. &
                     hy_conserveAngField) then
                     !! typically this should be zero. This is a hotfix...
                     !! do not use it with magnetic resistivity.
                     
                     
                   Sgeo(MAG_PHI) = 0.0!- ((pmomStar/densStar) * xmagStar - &
                                     !pmagStar * (xmomStar/densStar)) / xCenter(i) +&
                                     !(FR(MAG_PHI_FLUX)-FL(MAG_PHI_FLUX))/dx -&
                                     !(FR(MAG_PHI_FLUX)*((xCenter(i)+0.5*dx)/xCenter(i))-&
                                     ! FL(MAG_PHI_FLUX)*((xCenter(i)-0.5*dx)/xCenter(i)))/dx
                                     !if (rghtFac>0.0 .and. leftFac==0.0) print*,Sgeo(MAG_PHI)
                 endif
                 if (hy_geometry == CYLINDRICAL .and. &
                     hy_useMagneticResistivity) then
                   
                   !! This should not be a permanent solution! We need to
                   !! get the other induction formalism working (maybe by
                   !! splitting flux and press) so as to add the resisitive
                   !! correction to the FR FL dirrectly in addResistiveFluxes.
                   !! Now that flux is added at the RHSs which is plain ugly.
                   !! PT
                   
                   Sgeo(MAG_PHI) = Sgeo(MAG_PHI) + res_source(i,j,k) 
                 endif
#endif
                   Sgeo(MOM_PHI) = - Sgeo(MOM_PHI)


                   if (hy_geometry == SPHERICAL) then
                      tmomStar = U0(MOM_THT) 
                      Sgeo(HY_XMOM) = Sgeo(HY_XMOM) + tmomStar**2/densStar / xCenter(i)
                   endif
                endif

                call updateConservedVariable&
                     (U0(HY_DENS:HY_DENS+HY_VARINUM-1),&
                      FL(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1),&
                      FR(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1),&
                      GL(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1),&
                      GR(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1),&
                      HL(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1),&
                      HR(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1),&
                      gravX(i,j,k),gravY(i,j,k),gravZ(i,j,k),dx,dy,dz,dt,Sgeo,&
! <- ychen 10-2014
                      xcent(i), ycent(j), zcent(k))
! ychen ->

                U(DENS_VAR,i,j,k) = max(U0(HY_DENS),hy_smalldens)                    !density
#if defined(FLASH_USM_MHD) 
                U(MAGX_VAR:MAGZ_VAR,i,j,k) = U0(HY_MAGX:HY_MAGZ)                     !magnetic fields
#endif
                U(ENER_VAR,i,j,k) = U0(HY_ENER)                                      !total plasma energy
                !! We will update velocity fields after species & mass scalar update
#if (NSPECIES+NMASS_SCALARS) > 0
#ifdef FLASH_UHD_NEED_SCRATCHVARS
                !! Note that the velocity fields here are old velocities at time step n, not n+1
                call hy_uhd_updateSpeciesMassScalar&
                     (hy_order,U(DENS_VAR,i,j,k),&
                      SpOld(1:hy_numXN,i-3:i+3,j-3*ky:j+3*ky,k-3*kz:k+3*kz),&
                      Uold(1:6,        i-3:i+3,j-3*ky:j+3*ky,k-3*kz:k+3*kz),&
                      FL(HY_DENS_FLUX),FR(HY_DENS_FLUX),&
                      GL(HY_DENS_FLUX),GR(HY_DENS_FLUX),&
                      HL(HY_DENS_FLUX),HR(HY_DENS_FLUX),&
                      dx,dy,dz,dt, U(SPECIES_BEGIN:MASS_SCALARS_END,i,j,k))
#else
                sumSpecies = 0.
                do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
                   isph= ispu-NPROP_VARS
                   call updateSpeciesMassScalar&
                        (U(ispu,i,j,k),Uold(HY_DENS,i,j,k),U(DENS_VAR,i,j,k),&
                         FL(HY_END_FLUX+isph),FR(HY_END_FLUX+isph),&
                         GL(HY_END_FLUX+isph),GR(HY_END_FLUX+isph),&
                         HL(HY_END_FLUX+isph),HR(HY_END_FLUX+isph),&
                         dx,dy,dz,dt)

                   !! Conserving mass fractions
                   if (ispu <= SPECIES_END) then
                      sumSpecies = sumSpecies + U(ispu,i,j,k)
                   endif

                enddo

                !! Conserving mass fractions
                do ispu = SPECIES_BEGIN, SPECIES_END
                   isph = ispu-NPROP_VARS
                   if (ispu <= SPECIES_END) then
                      U(ispu,i,j,k) = U(ispu,i,j,k)/sumSpecies
                   endif
                enddo

#endif /*ifdef FLASH_UHD_NEED_SCRATCHVARS */
#endif /*if (NSPECIES+NMASS_SCALARS) > 0*/

                !! Update velocity fields after species & mass scalar update
                U(VELX_VAR:VELZ_VAR,i,j,k) = U0(HY_XMOM:HY_ZMOM)/U(DENS_VAR,i,j,k)        !velocities

#ifdef FLASH_UHD_3T
#ifdef FLASH_UHD_HYDRO
                !! Perform energy updates here for 3T.  For 1T energy is updated in hy_uhd_energyFix.
                !! We only perform the energy updates here for pure hydro only.
                !! MHD still needs to update magnetic fields after the current cell-centered
                !! variable udpates. Therefore, the energy updates should be done after
                !! updating magnetic fields.
                ekin = .5*dot_product(U(VELX_VAR:VELZ_VAR,i,j,k),U(VELX_VAR:VELZ_VAR,i,j,k))&
                     *U(DENS_VAR,i,j,k)

                eint = U(ENER_VAR,i,j,k)-ekin

                if (.not. hy_useAuxEintEqn .or. eint > hy_eswitch*ekin) then
                   U(EINT_VAR,i,j,k) = max(hy_smallE,eint/U(DENS_VAR,i,j,k))
                else
                   U(EINT_VAR,i,j,k) = max(hy_smallE,U(EINT_VAR,i,j,k)/U(DENS_VAR,i,j,k))
                endif
                !! Store specific gas energy ener = ekin + eint
                U(ENER_VAR,i,j,k) = U(EINT_VAR,i,j,k) + ekin/U(DENS_VAR,i,j,k)
#endif /*ifdef FLASH_UHD_HYDRO*/
#endif /*ifdef FLASH_UHD_3T*/

#ifdef BDRY_VAR
             endif
#endif

          enddo !end of i loop
       enddo !end of j loop
    enddo !end of k loop

    !! ---------------------------------------------------------------
#if defined(FLASH_USM_MHD) 
    if (hy_forceHydroLimit) then
       U(MAGX_VAR:MAGZ_VAR,:,:,:) = 0.
    endif
#endif

    ! Release block pointers
    call Grid_releaseBlkPtr(blockID,U,CENTER)
    call Grid_releaseBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)

  End Subroutine hy_uhd_unsplitUpdate



!! ==================================================================
Subroutine updateConservedVariable(Ul,FL,FR,GL,GR,HL,HR,  &
                                      gravX, gravY, gravZ,&
                                      dx,dy,dz,dt,Sgeo, &
! <- ychen 10-2014
                                      x,y,z)
! ychen ->
  use Hydro_data, ONLY : hy_useGravity
! <- ychen 10-2014
  use Simulation_data
  use Driver_data,  ONLY : dr_simTime
! ychen ->
  implicit none
  real, dimension(HY_VARINUM), intent(INOUT) :: Ul,Sgeo
  real, dimension(HY_VARINUM), intent(IN) :: FL,FR,GL,GR,HL,HR
  real, intent(IN) :: gravX,gravY,gravZ
  real, intent(IN) :: dx,dy,dz,dt
  real, dimension(3) :: momentaOld
  real :: densOld
! <- ychen 10-2014
  real, dimension(HY_VARINUM) :: UdeltaNozzle
  real, intent(IN) :: x, y, z
  integer :: nozzle=1
  real :: radius, length, sig, distance, theta, vel, fac
  real, dimension(3) :: cellvec, plnvec, jetvec, rvec, phivec, velvec

  cellvec = (/ x, y, z /)
  call hy_uhd_jetNozzleGeometry(nozzle,cellvec,radius,length,distance,&
                                sig,theta,jetvec,rvec,plnvec,phivec)
! ychen ->

  !! Store old states at n
  densOld = Ul(HY_DENS)
  momentaOld(1:3) = Ul(HY_XMOM:HY_ZMOM)

! <- ychen 10-2014
  if ((radius.le.(sim(nozzle)%radius+sim(nozzle)%rfeather_outer))&
      .and.(abs(length).le.(sim(nozzle)%length+sim(nozzle)%zfeather))) then
     ! inside the jet nozzle
     !! Update conservative variables from n to n+1 step
     vel = sim(nozzle)%velocity*sin(PI/2.0*min(abs(length),sim(nozzle)%length)*sig/sim(nozzle)%length)
     fac = taper(nozzle, radius, length, 1.0, 1.0, 0.0)
     velvec = vel*jetvec + 0.1*sim(nozzle)%velocity*plnvec*&
                      0.5*(1.0+cos(PI*(min(0.0, radius-sim(nozzle)%radius)/sim(nozzle)%rfeather_outer)))
     !fac = taperR(nozzle, radius, 1.0, 0.0)
     ! smooth transition from nozzle to flash grid
     ! 1.0 for nozzle injection; 0.0 for flash solution
     UdeltaNozzle(HY_DENS)        =sim(nozzle)%deltaRho*fac
     UdeltaNozzle(HY_XMOM:HY_ZMOM)=velvec*fac*sim(nozzle)%deltaRho*fac
     UdeltaNozzle(HY_ENER)        =sim(nozzle)%deltaP*fac/(sim(nozzle)%gamma-1.0)+&
                                   0.5*sim(nozzle)%deltaRho*fac*(vel*fac)**2
     UdeltaNozzle(HY_MAGX:HY_MAGZ)=(/0.0, 0.0, 0.0/)
     !write (*,'(a10, 10e9.2, f6.3)') 'Unozzle', Unozzle, radius, length, fac
     !write (*,'(a10, 8e9.2)') 'Ul', Ul

     Ul(HY_DENS:HY_DENS+HY_VARINUM-1)=Ul(HY_DENS:HY_DENS+HY_VARINUM-1)&
          +UdeltaNozzle(HY_DENS:HY_DENS+HY_VARINUM-1)&
          -dt/dx*( FR(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1)&
                  -FL(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1))*(1.0-fac)
     !write (*,'(a10, 8e9.2)') 'Ulx', Ul
     
     if (NDIM > 1) then
     Ul(HY_DENS:HY_DENS+HY_VARINUM-1)=Ul(HY_DENS:HY_DENS+HY_VARINUM-1)&
          -dt/dy*( GR(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1)&
                  -GL(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1))*(1.0-fac)
     !write (*,'(a10, 8e9.2)') 'Uly', Ul

        if (NDIM > 2) then
        Ul(HY_DENS:HY_DENS+HY_VARINUM-1)=Ul(HY_DENS:HY_DENS+HY_VARINUM-1)& 
             -dt/dz*( HR(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1)&
                     -HL(HY_DENS_FLUX:HY_DENS_FLUX+HY_VARINUM-1))*(1.0-fac)
        endif

     endif
     !write (*,'(a10, 8e9.2)') 'HR', HR
     !write (*,'(a10, 8e9.2)') 'HL', HL
     !write (*,'(a10, 8e9.2)') 'Ulz', Ul
  else
! ychen ->
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
! <- ychen 10-2014
  endif
! ychen ->


  !! Include geometric source term
  Ul = Ul + dt*Sgeo

  if (hy_useGravity) then
     !! Note: 1. Extrapolated gravity case
     !!          gravX,Y,Z - gravity at n
     !!------------------------------------------

     !! The following is a new gravity formulation:
     !! Here we only update total energy (rho*E) with gravity source at n.
     !! The final gravity update is to be done in hy_uhd_unsplit 
     !! after this conservative upate this routine.
     !! The second part includes gravity at a new state n+1, finalizing the
     !! gravity coupling to hydro/MHD variables.
     Ul(HY_XMOM:HY_ZMOM) = Ul(HY_XMOM:HY_ZMOM)&
          + 0.5*dt*densOld*(/gravX,gravY,gravZ/)

     Ul(HY_ENER) = Ul(HY_ENER) &
          + 0.5*dt*dot_product(momentaOld(1:3),(/gravX,gravY,gravZ/))

  endif


End Subroutine updateConservedVariable


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
