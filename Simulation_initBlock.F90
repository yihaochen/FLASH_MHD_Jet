!!****if* source/Simulation/SimulationMain/WindTunnel/Simulation_initBlock
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
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Mach 3 wind tunnel
!!  problem.
!!
!!  References:  Emery, A. E., 1968, JCP, 2, 306
!!               Woodward, P. and Colella, P., 1984, JCP, 54, 115
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
  use Logfile_interface, ONLY : Logfile_stamp
  use Simulation_data
  use Driver_data, ONLY : dr_simTime


  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) :: blockID
  real,pointer,dimension(:,:,:,:) :: solnData, &
        solnFaceXData, solnFaceYData, solnFaceZData


  real :: rho_zone, velx_zone, vely_zone, velz_zone, pres_zone, &
       ener_zone, ekin_zone, eint_zone, vel


  integer :: i,j,k, istat
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  logical :: gcell = .true.

  logical, save :: firstCall = .TRUE.

  integer :: nozzle=1
  real, dimension(3) :: cellvec
  real :: radius, length, sig, distance, theta
  real, dimension(3) :: plnvec, jetvec, rvec, phivec

!===============================================================================

  if (firstCall) then


  endif

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
  

! In this problem the initial conditions are spatially uniform.
  
  rho_zone = sim_rhoAmbient
  pres_zone = sim_pAmbient

  velx_zone = sim_windVel
  vely_zone = 0.0
  velz_zone = 0.0


  ! Compute the gas energy and set the gamma-values needed for
  ! the equation of state.
  ekin_zone = 0.5 * (velx_zone**2 + vely_zone**2 + velz_zone**2)

  eint_zone = pres_zone / (sim_gamma-1.) / rho_zone
  ener_zone = eint_zone + ekin_zone
  ener_zone = max(ener_zone, sim_smallP)


#if NSPECIES > 0
  solnData(SPECIES_BEGIN,:,:,:) =  1.0-(NSPECIES-1)*sim_smallX
  solnData(SPECIES_BEGIN+1:SPECIES_END,:,:,:) =     sim_smallX
#endif

  ! store the variables in the block's unk data
  solnData(DENS_VAR,:,:,:) = rho_zone
  solnData(PRES_VAR,:,:,:) = pres_zone
  solnData(ENER_VAR,:,:,:) = ener_zone
#ifdef EINT_VAR
  solnData(EINT_VAR,:,:,:) = eint_zone
#endif
  solnData(GAMC_VAR,:,:,:) = sim_gamma
  solnData(GAME_VAR,:,:,:) = sim_gamma


  solnData(VELX_VAR,:,:,:) = velx_zone
  solnData(VELY_VAR,:,:,:) = vely_zone
  solnData(VELZ_VAR,:,:,:) = velz_zone

  solndata(MAGX_VAR,:,:,:) = 0.0
  solndata(MAGY_VAR,:,:,:) = 0.0
  solndata(MAGZ_VAR,:,:,:) = sim_bzAmbient
  solndata(MAGP_VAR,:,:,:) = max(sim_bzAmbient*sim_bzAmbient/8.0/PI, sim_smallP)

  solnFaceXdata(MAG_FACE_VAR,:,:,:) = 0.0
  solnFaceYdata(MAG_FACE_VAR,:,:,:) = 0.0
  solnFaceZdata(MAG_FACE_VAR,:,:,:) = sim_bzAmbient

  ! Initialize the nozzle
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
   do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
       cellvec = (/ xCoord(i), yCoord(j), zCoord(k) /)
       call hy_uhd_jetNozzleGeometry(nozzle,cellvec,radius,length,distance,sig,theta,jetvec,rvec,plnvec,phivec)
       if ((radius.le.2.0*(sim(nozzle)%radius+sim(nozzle)%bfeather_outer))&
           .and.(abs(length).le.2.0*(sim(nozzle)%length+sim(nozzle)%zfeather))) then
          !if ((radius.gt.2.0*sim(nozzle)%radius).or.(abs(length).gt.2.0*sim(nozzle)%length)) then
             solnData(DENS_VAR,i,j,k) = taper(nozzle, 0.5*radius, 0.5*length, sim(nozzle)%density,&
                                              sim(nozzle)%density, sim_rhoAmbient)
             solnData(PRES_VAR,i,j,k) = taper(nozzle, 0.5*radius, 0.5*length, sim(nozzle)%pressure,&
                                              sim(nozzle)%pressure, sim_pAmbient)
             if (abs(length).gt.(sim(nozzle)%length)) then
                 vel = taper(nozzle, radius, 0.5*length, sim(nozzle)%velocity,&
                             sim(nozzle)%velocity, 0.0)*sig
             else
                 vel = taper(nozzle, radius, 0.5*length, sim(nozzle)%velocity,&
                             sim(nozzle)%velocity, 0.0)*sin(PI/2.0*length/sim(nozzle)%length)
             endif
             solnData(VELX_VAR,i,j,k) = vel*jetvec(1)
             solnData(VELY_VAR,i,j,k) = vel*jetvec(2)
             solnData(VELZ_VAR,i,j,k) = vel*jetvec(3)
             solnData(JET_SPEC,i,j,k) = taper(nozzle, 0.5*radius, 0.5*length, 1.0, 1.0, 0.0)
             solnData(ISM_SPEC,i,j,k) = taper(nozzle, 0.5*radius, 0.5*length, 0.0, 0.0, 1.0)
          !else
          !   vel = sim(nozzle)%velocity*sin(PI/2.0*length/sim(nozzle)%length)
          !   solnData(DENS_VAR,i,j,k) = sim(nozzle)%density
          !   solnData(PRES_VAR,i,j,k) = sim(nozzle)%pressure
          !   solnData(VELX_VAR,i,j,k) = vel*jetvec(1)*sig
          !   solnData(VELY_VAR,i,j,k) = vel*jetvec(2)*sig
          !   solnData(VELZ_VAR,i,j,k) = vel*jetvec(3)*sig
          !   solnData(JET_SPEC,i,j,k) = 1.0 - sim_smallX
          !   solnData(ISM_SPEC,i,j,k) = sim_smallX
          !endif
          solnData(EINT_VAR,i,j,k) = solnData(PRES_VAR,i,j,k)/solnData(DENS_VAR,i,j,k)&
                                     /(solnData(GAME_VAR,i,j,k)-1.0)
          solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)+&
                                0.5*(solnData(VELX_VAR,i,j,k)**2 +&
                                     solnData(VELY_VAR,i,j,k)**2 +&
                                     solnData(VELZ_VAR,i,j,k)**2)
       endif
    enddo
   enddo
  enddo
  call Heat_fillnozzle(blockID, 0.0, dr_simTime, .true.)

  call Grid_releaseBlkPtr(blockID, solnData, CENTER)
  call Grid_releaseBlkPtr(blockID,solnFaceXData,FACEX)
  call Grid_releaseBlkPtr(blockID,solnFaceYData,FACEY)
  call Grid_releaseBlkPtr(blockID,solnFaceZData,FACEZ)
 

  return
end subroutine Simulation_initBlock



