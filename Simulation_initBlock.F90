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
       ener_zone, ekin_zone, eint_zone


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


  call Grid_getBlkPtr(blockID,solnData,CENTER)
  call Grid_getBlkPtr(blockID,solnFaceXData,FACEX)
  call Grid_getBlkPtr(blockID,solnFaceYData,FACEY)
  call Grid_getBlkPtr(blockID,solnFaceZData,FACEZ)

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
  solndata(MAGZ_VAR,:,:,:) = 0.0
  solndata(MAGP_VAR,:,:,:) = sim_smallP

  solnFaceXdata(MAG_FACE_VAR,:,:,:) = 0.0
  solnFaceYdata(MAG_FACE_VAR,:,:,:) = 0.0
  solnFaceZdata(MAG_FACE_VAR,:,:,:) = 0.0

  ! Initialize the nozzle
  call Heat_fillnozzle(blockID, 0.0, dr_simTime)

  call Grid_releaseBlkPtr(blockID, solnData, CENTER)
  call Grid_releaseBlkPtr(blockID,solnFaceXData,FACEX)
  call Grid_releaseBlkPtr(blockID,solnFaceYData,FACEY)
  call Grid_releaseBlkPtr(blockID,solnFaceZData,FACEZ)
 

  return
end subroutine Simulation_initBlock



