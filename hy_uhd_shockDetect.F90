!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_shockDetect
!!
!! NAME
!!
!!  hy_uhd_shockDetect
!!
!! SYNOPSIS
!!
!!  hy_uhd_shockDetect( integer (IN) :: blockID )
!!
!! DESCRIPTION
!!
!!  This routine detects strongly compressive motions in simulation
!!  by calculating undivided pressure gradients and divergence of
!!  velocity fields. Two parameters beta and delta have been set 
!!  to detect strong shocks. If such shocks exist then the unsplit
!!  scheme applies its robust flux differencings using limited slopes
!!  in data reconstruction step (see hy_uhd_dataReconstruct.F90).
!!  Different shock strengths can also be detected by lowering/increasing
!!  beta and delta values.
!!
!! ARGUMENTS
!!
!!  blockID  - local block ID
!!
!! REFERENCE 
!!
!!  Balsara and Spicer, JCP, 149:270--292, 1999.
!!
!!***

!!REORDER(4): U

Subroutine hy_uhd_shockDetect( blockID )


  use Grid_interface,    ONLY : Grid_getBlkIndexLimits, &
                                Grid_getBlkPtr,         &
                                Grid_releaseBlkPtr,     &
                                Grid_getDeltas
  use Hydro_data,        ONLY : hy_cfl, hy_cfl_original,&
                                hy_RiemannSolver,       &
                                hy_geometry,            &
                                hy_fallbackLowerCFL,    &
                                hy_shockLowerCFL
  use Driver_data,       ONLY : dr_nStep
  use Logfile_interface, ONLY : Logfile_open,Logfile_close
  use Driver_interface,  ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockID
  !! -----------------------------------------------------

  integer :: i,j,k
  logical :: SW1, SW2
  logical :: doLowerCFL
  integer :: k2,k3

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC

  real :: divv,gradPx,gradPy,gradPz
  real :: minP,minC,beta,delta
  real :: localCfl,cflMax
  real, dimension(:,:,:), allocatable :: Vc
  real, dimension(:,:,:,:), pointer   :: U


  ! Two parameters that can be adjusted to detect shocks
  ! with different strengths:
  ! (a) The lower the values the weaker shocks detected 
  !     (lower the values to detect more shock regions)
  ! (b) The larger the values the stronger shocks detected
  !     (increase the values to detect less shock regions)
  beta = 0.5 !0.5 !10. ! gradP
  delta= 0.1  !0.1 !2. ! divV
!!$  beta  = 0.1 !0.1
!!$  delta = 0.01

  k2=0
  k3=0
  ! Set dimensional indices
  if (NDIM > 1) k2=1
  if (NDIM > 2) k3=1


  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getBlkPtr(blockID,U,CENTER)

#ifdef SHOK_VAR
     U(SHOK_VAR,:,:,:)=0.
#else
  if (hy_RiemannSolver == HYBR) then
     call Driver_abortFlash&
          ("[hy_uhd_shockDetect]: SHOK_VAR has not been defined for shock detection")
  endif
#endif


  !! Allocate a temporary cell-centered array for sound speed
  allocate(Vc(blkLimits(LOW,IAXIS)- 1:blkLimits(HIGH,IAXIS)+1,  &
              blkLimits(LOW,JAXIS)-k2:blkLimits(HIGH,JAXIS)+k2, &
              blkLimits(LOW,KAXIS)-k3:blkLimits(HIGH,KAXIS)+k3))

  !! Compute sound speed
  do k=blkLimits(LOW,KAXIS)-K3D,blkLimits(HIGH,KAXIS)+K3D
     do j=blkLimits(LOW,JAXIS)-K2D,blkLimits(HIGH,JAXIS)+K2D
        do i=blkLimits(LOW,IAXIS)-1,blkLimits(HIGH,IAXIS)+1
           Vc(i,j,k) = sqrt(U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k))
        enddo
     enddo
  enddo

  ! Note: it is the job of the caller to revert back to the original cfl, if desired
  doLowerCFL = .FALSE.

  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           localCfl = hy_cfl

           ! initialize switch values
           SW1 = .false.
           SW2 = .false.

           minP = minval(U(PRES_VAR,i-1:i+1,j-k2:j+k2,k-k3:k+k3))
           minC = minval(        Vc(i-1:i+1,j-k2:j+k2,k-k3:k+k3))

           !! We do not need to include non-Cartesian geom factors here.
           !! Undivided divV
           divv =        U(VELX_VAR,i+1,j,  k  ) - U(VELX_VAR,i-1,j,  k  )
#if NDIM > 1
           divv = divv + U(VELY_VAR,i,  j+1,k  ) - U(VELY_VAR,i,  j-1,k  )
#if NDIM == 3
           divv = divv + U(VELZ_VAR,i,  j,  k+1) - U(VELZ_VAR,i,  j,  k-1)
#endif
#endif
           divv = 0.5*divv  

           !! Undivided grad pres
           gradPx = 0.5*(U(PRES_VAR,i+1,j,  k  ) - U(PRES_VAR,i-1,j,  k  ))
           gradPy = 0.
           gradPz = 0.
#if NDIM > 1
           gradPy = 0.5*(U(PRES_VAR,i,  j+1,k  ) - U(PRES_VAR,i,  j-1,k  ))
#if NDIM == 3
           gradPz = 0.5*(U(PRES_VAR,i,  j,  k+1) - U(PRES_VAR,i,  j,  k-1))
#endif
#endif
           if ( abs(gradPx)+abs(gradPy)+abs(gradPz) .ge. beta*minP ) then
              SW1 = .true.
           endif

           if (-delta*minC .ge. divv) then
              SW2 = .true.
           endif


           if (SW1 .and. SW2) then
              ! Set SHOCK_VAR to 1.0 if a shock is detected. 
              ! One use is for a local hybrid method in the Hydro unit which
              ! applies (a diffusive) HLL solver when SHOK_VAR = 1.
#ifdef SHOK_VAR
              U(SHOK_VAR,i,j,k) = 1.
#endif

              if (hy_shockLowerCFL) then
                 ! if lowering of the CFL factor within shocks is requested...

                 doLowerCFL = .TRUE.

#ifdef CFL_VAR
                 if (NDIM == 1) then
                    if (localCfl > 0.60) localCfl = 0.60
                 elseif (NDIM == 2) then
                    if (localCfl > 0.45) localCfl = 0.45
                 elseif (NDIM == 3) then
                    if (localCfl > 0.25) localCfl = 0.25
                 endif
#endif
              endif ! endif (hy_shockLowerCFL) then

           endif !endif (SW1 .and. SW2) then

#ifdef CFL_VAR
           if (hy_shockLowerCFL .OR. .NOT. hy_fallbackLowerCFL) &
                U(CFL_VAR,i,j,k) = localCfl
#endif
        enddo !enddo i-loop
     enddo !enddo j-loop
  enddo !enddo k-loop

  ! Release block pointer
  call Grid_releaseBlkPtr(blockID,U,CENTER)

  ! Deallocate sound speed array
  deallocate(Vc)

  if (doLowerCFL) then
     if (NDIM == 1) then
        cflMax = 0.60
     elseif (NDIM == 2) then
        cflMax = 0.45
     elseif (NDIM == 3) then
        cflMax = 0.25
     endif

!!$     !$omp critical (Update_cfl)
!!$     hy_cfl = min(hy_cfl,cflMax)
!!$     !$omp end critical (Update_cfl)
!!$
     !$omp atomic 
     hy_cfl = min(hy_cfl,cflMax)

  endif

End Subroutine hy_uhd_shockDetect
