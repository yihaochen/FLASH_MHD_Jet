!!****if* source/Simulation/SimulationMain/WindTunnel/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID, 
!!                       integer(IN) :: myPE)
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
!!  myPE -             my processor number
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(blockID, myPE)

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getCellCoords
  use Logfile_interface, ONLY : Logfile_stamp
  use Simulation_data


  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) :: blockID, myPE
  real,pointer :: solnData(:,:,:,:)

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

     if (myPE == MASTER_PE) then

        call Logfile_stamp(myPE,"initializing for jet ", 'run_init')
        write (*,*) "flash:  initializing for jet"

     endif


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

  eint_zone = pres_zone / (sim_gamma-1.)
  eint_zone = eint_zone / rho_zone
  ener_zone = eint_zone + ekin_zone
  ener_zone = max(ener_zone, sim_smallP)


  call Grid_getBlkPtr(blockID, solnData, CENTER)
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

  ! Initialize the nozzle
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX),stat=istat)
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY),stat=istat)
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ),stat=istat)
  
  call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell, zCoord, sizeZ)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell, xCoord, sizeX)
  
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH, IAXIS)
           cellvec = (/ xCoord(i), yCoord(j), zCoord(k) /)
           call hy_uhd_jetNozzleGeometry(nozzle,cellvec,radius,length,distance,sig,theta,jetvec,rvec,plnvec,phivec)
           if ( radius.lt.sim(nozzle)%radius.and.(abs(length).lt.sim(nozzle)%length) ) then
           ! inside the jet nozzle
              solnData(VELX_VAR,i,j,k) = sim(nozzle)%velocity*jetvec(1)*sig
              solnData(VELY_VAR,i,j,k) = sim(nozzle)%velocity*jetvec(2)*sig
              solnData(VELZ_VAR,i,j,k) = sim(nozzle)%velocity*jetvec(3)*sig
              solnData(DENS_VAR,i,j,k) = sim(nozzle)%density
              solnData(PRES_VAR,i,j,k) = sim(nozzle)%pressure
              solnData(EINT_VAR,i,j,k) = 1.0/(solnData(GAME_VAR,i,j,k)-1.0)
              solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)+&
                                    0.5*(solnData(VELX_VAR,i,j,k)**2 + &
                                         solnData(VELY_VAR,i,j,k)**2 + &
                                         solnData(VELZ_VAR,i,j,k)**2)
              solnData(JET_SPEC,i,j,k) = 1.0 - sim_smallX
              solnData(ISM_SPEC,i,j,k) = sim_smallX
           endif
        enddo
     enddo
  enddo

  ! force refinement to ensure the nozzle is resolved
  ! not robust, shuold be replaced by gr_markWithRadius.F90
  !if (firstCall) then
  !   solnData(DENS_VAR,1,1,1) = sim(nozzle)%density
  !endif
  firstCall = .FALSE.

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

  call Grid_releaseBlkPtr(blockID, solnData, CENTER)
 

  return
end subroutine Simulation_initBlock



