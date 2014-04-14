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
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
    Grid_getCellCoords, Grid_putPointData
  !use Hydro_data, ONLY: hy_unsplitEosMode
  !use Eos_interface, ONLY : Eos_wrapped
  use Simulation_data
  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN)::blockList
  real,intent(IN) :: dt,time

  integer :: blockID, i,j,k, istat
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  logical :: gcell = .true.
  real, pointer, dimension(:,:,:,:) :: solnData

  integer :: nozzle=1
  real, dimension(3) :: cellvec
  real :: radius, length, sig, distance, theta
  real, dimension(3) :: plnvec, jetvec, rvec, phivec

  do blockID=1,blockCount
     
     call Grid_getBlkIndexLimits(blockList(blockID),blkLimits,blkLimitsGC)
     sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
     allocate(xCoord(sizeX),stat=istat)
     sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
     allocate(yCoord(sizeY),stat=istat)
     sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
     allocate(zCoord(sizeZ),stat=istat)
  
     call Grid_getCellCoords(KAXIS,blockList(blockID),CENTER,gcell, zCoord, sizeZ)
     call Grid_getCellCoords(JAXIS,blockList(blockID),CENTER,gcell, yCoord, sizeY)
     call Grid_getCellCoords(IAXIS,blockList(blockID),CENTER,gcell, xCoord, sizeX)
     
     call Grid_getBlkPtr(blockList(blockID),solnData,CENTER)

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
              elseif ( (sim_accx.gt.sim_smallX).or.(sim_accy.gt.sim_smallX) ) then
              ! apply gravity
                 solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) + dt*sim_accx
                 solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) + dt*sim_accy
                 solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)+&
                                       0.5*(solnData(VELX_VAR,i,j,k)**2 + &
                                            solnData(VELY_VAR,i,j,k)**2 + &
                                            solnData(VELZ_VAR,i,j,k)**2)
              endif
           enddo
        enddo
     enddo

     deallocate(xCoord)
     deallocate(yCoord)
     deallocate(zCoord)

     !call Eos_wrapped(hy_unsplitEosMode, blkLimits, blockList(blockID))

     call Grid_releaseBlkPtr(blockList(blockID),solnData,CENTER)

  enddo

  return

end subroutine Heat


