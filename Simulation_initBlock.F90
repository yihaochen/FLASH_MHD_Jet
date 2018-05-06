!!****if* source/Simulation/SimulationMain/magnetoHD/MHD_Jet/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID)
!!
!!
!!
!! DESCRIPTION
!!
!!
!!
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(blockID)

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getCellCoords
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Simulation_data
  use Driver_data, ONLY : dr_simTime, dr_dtInit
  use Hydro_data, ONLY : hy_bref


  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) :: blockID
  real,pointer,dimension(:,:,:,:) :: solnData, &
        solnFaceXData, solnFaceYData, solnFaceZData


  real :: rho_zone, velx_zone, vely_zone, velz_zone, temp_zone, &
       ener_zone, ekin_zone, eint_zone, gasConst
  real :: densityBG, tempBG!, rhoCut
  real :: vel, fac


  integer :: i,j,k, istat
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  logical :: gcell = .true.

  integer :: nozzle=1
  real :: radius, length, sig, distance, theta
  real, dimension(3) :: cellvec, plnvec, jetvec, rvec, phivec, velvec, voutvec
  real :: r2, bf, rout, rmix

!===============================================================================

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(xCoord(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell, xCoord, sizeX)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell, zCoord, sizeZ)

  call Grid_getBlkPtr(blockID,solnData,CENTER)
  call Grid_getBlkPtr(blockID,solnFaceXData,FACEX)
  call Grid_getBlkPtr(blockID,solnFaceYData,FACEY)
  call Grid_getBlkPtr(blockID,solnFaceZData,FACEZ)
  

  ! Initial conditions are spatially uniform.

  call PhysicalConstants_get("ideal gas constant", gasConst)
  
  rho_zone = sim_rhoCore
  temp_zone = sim_Tcore

  velx_zone = 0.0
  vely_zone = sim_windVel
  velz_zone = 0.0


  ! Compute the gas energy and set the gamma-values needed for
  ! the equation of state.
  !ekin_zone = 0.5 * (velx_zone**2 + vely_zone**2 + velz_zone**2)

  !eint_zone = pres_zone / (sim_gamma-1.) / rho_zone
  !ener_zone = eint_zone + ekin_zone
  !ener_zone = max(ener_zone, sim_smallP)


#if NSPECIES > 0
  solnData(SPECIES_BEGIN,:,:,:) =  1.0-(NSPECIES-1)*sim_smallx
  solnData(SPECIES_BEGIN+1:SPECIES_END,:,:,:) =     sim_smallx
#endif

  ! store the variables in the block's unk data
  solnData(DENS_VAR,:,:,:) = rho_zone
  solnData(TEMP_VAR,:,:,:) = temp_zone
  solnData(GAMC_VAR,:,:,:) = sim_gamma
  solnData(GAMC_VAR,:,:,:) = sim_gamma
  solnData(VELX_VAR,:,:,:) = velx_zone
  solnData(VELY_VAR,:,:,:) = vely_zone
  solnData(VELZ_VAR,:,:,:) = velz_zone

#ifdef EINT_VAR
  !solnData(EINT_VAR,:,:,:) = eint_zone
  !solnData(EINT_VAR,:,:,:) = solnData(PRES_VAR,:,:,:)/solnData(DENS_VAR,:,:,:)&
  !                                   /(solnData(GAMC_VAR,:,:,:)-1.0)
  solnData(EINT_VAR,:,:,:) = gasConst*solnData(TEMP_VAR,:,:,:)&
                             /(solnData(GAMC_VAR,:,:,:)-1.0)/sim_mu
#endif
  !solnData(ENER_VAR,:,:,:) = ener_zone
  solnData(ENER_VAR,:,:,:) = solnData(EINT_VAR,:,:,:)+&
                        0.5*(solnData(VELX_VAR,:,:,:)**2 +&
                             solnData(VELY_VAR,:,:,:)**2 +&
                             solnData(VELZ_VAR,:,:,:)**2)

  solndata(MAGX_VAR,:,:,:) = 0.0
  solndata(MAGY_VAR,:,:,:) = 0.0
  solndata(MAGZ_VAR,:,:,:) = sim_bzAmbient
  solndata(MAGP_VAR,:,:,:) = max((sim_bzAmbient/hy_bref)**2/2.0, sim_smallP)

  solnFaceXdata(MAG_FACE_VAR,:,:,:) = 0.0
  solnFaceYdata(MAG_FACE_VAR,:,:,:) = 0.0
  solnFaceZdata(MAG_FACE_VAR,:,:,:) = sim_bzAmbient

  !rhoCut = sim_rhoCore*(1.0 + (sim_rCut/sim_rCore)**2)**(-1.5*sim_densityBeta)

  bf = sim(nozzle)%rFeatherOut
  r2 = sim(nozzle)%radius
  rout = r2 + bf
  rmix = rout + sim(nozzle)%rFeatherMix

  ! Initialize the nozzle and environment
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
   do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
    do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
       cellvec = (/ xCoord(i), yCoord(j), zCoord(k) /)
       call hy_uhd_jetNozzleGeometry(nozzle,cellvec,radius,length,distance,&
                                     sig,theta,jetvec,rvec,plnvec,phivec)
       ! inside the nozzle
       if ((radius.le.rmix)&
           .and.(abs(length).le.2.0*(sim(nozzle)%length+sim(nozzle)%zFeather))) then
          fac = taper(nozzle, radius, 0.5*length, 1.0, 1.0, 0.0)
          !vel = sim(nozzle)%velocity&
          !      *(0.5*(1.0+cos(PI*(max(0.0, min(1.0, (radius-r2)/bf)))))&
          !      *(1.0-sim(nozzle)%outflowR)+sim(nozzle)%outflowR ) &
          !      *sin(PI/2.0*min(abs(length),0.5*sim(nozzle)%length)*sig/sim(nozzle)%length/0.5)
          !voutvec = sim(nozzle)%outflowR*sim(nozzle)%velocity*plnvec&
          !          !*coshat(radius-0.5*(r2+2.0*bf), 0.5*(r2+bf), bf, 1.0)
          !          *0.5*(1.0+cos(PI*( max(-1.0, min(0.0,(radius-rout)/bf)) )))

          !velvec = vel*jetvec + voutvec &
          !         + sim(nozzle)%linVel + cross(sim(nozzle)%angVel,rvec*distance)
          !solnData(VELX_VAR:VELZ_VAR,i,j,k) = velvec*fac + solnData(VELX_VAR:VELZ_VAR,i,j,k)*(1.0-fac)
       endif
       ! cylindrical initial cavity
       if (sim(nozzle)%initGeometry == 'cylindrical') then
          if ((radius.le.rmix)&
              .and.(abs(length).le.2.0*(sim(nozzle)%length+sim(nozzle)%zFeather))) then
             ! inside the extended nozzle and feather
             fac = taper(nozzle, radius, 0.5*length, 1.0, 1.0, 0.0)
             !fac = 1.0
          else
             fac = 0.0
          endif
       ! spherical initial cavity
       else
          if (distance.le.2.0*max( rmix, (sim(nozzle)%length+sim(nozzle)%zFeather) ) ) then
             fac = taperSph(nozzle, 0.5*distance, 1.0, 0.0)
          else
             fac = 0.0
          endif
       endif

       ! Background
       if (sim_densityProfile =="betacore") then
          ! beta model for the density profile
          densityBG = max(sim_rhoCore*(1.0 + (distance/sim_rCore)**2)**(-1.5*sim_densityBeta),&
                          sim_rhoFloor)
          ! isothermal two-temperature atmosphere
          tempBG = sim_Tout*(1.0+(distance/sim_rCoreT)**3)&
                       /(sim_Tout/sim_Tcore+(distance/sim_rCoreT)**3)
       else
          ! uniform background density and pressure
          densityBG = sim_rhoCore
          tempBG = sim_Tcore
       endif
       solnData(DENS_VAR,i,j,k) = sim(nozzle)%density*fac + densityBG*(1.0-fac)
       solnData(PRES_VAR,i,j,k) = sim(nozzle)%pressure*fac &
                                  + tempBG*densityBG*gasConst/sim_mu*(1.0-fac)
       !solnData(TEMP_VAR,i,j,k) = sim(nozzle)%pressure/sim(nozzle)%density/gasConst*sim_mu*fac &
       !                           + tempBG*(1.0-fac)
       solnData(JET_SPEC,i,j,k) = max(sim_smallx, fac)
       solnData(ISM_SPEC,i,j,k) = max(sim_smallx, 1.0-fac)
       !solnData(EINT_VAR,i,j,k) = max(sim_smalle, solnData(PRES_VAR,i,j,k)/solnData(DENS_VAR,i,j,k)&
       !                           /(solnData(GAMC_VAR,i,j,k)-1.0))
       !solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)+&
       !                      0.5*(solnData(VELX_VAR,i,j,k)**2 +&
       !                           solnData(VELY_VAR,i,j,k)**2 +&
       !                           solnData(VELZ_VAR,i,j,k)**2)
    enddo
   enddo
  enddo

  call Eos_wrapped(MODE_DENS_PRES, blkLimitsGC, blockID, CENTER)
  solnData(ENER_VAR,:,:,:) = solnData(EINT_VAR,:,:,:)+&
                        0.5*(solnData(VELX_VAR,:,:,:)**2 +&
                             solnData(VELY_VAR,:,:,:)**2 +&
                             solnData(VELZ_VAR,:,:,:)**2)
  call Grid_releaseBlkPtr(blockID, solnData, CENTER)
  call Grid_releaseBlkPtr(blockID,solnFaceXData,FACEX)
  call Grid_releaseBlkPtr(blockID,solnFaceYData,FACEY)
  call Grid_releaseBlkPtr(blockID,solnFaceZData,FACEZ)
 

  return
end subroutine Simulation_initBlock



