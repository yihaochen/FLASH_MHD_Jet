!!****f* source/physics/sourceTerms/Heat/Heat
!!
!! NAME
!!  
!!  Heat_fillnozzle 
!!
!!
!! SYNOPSIS
!! 
!!  call Heat_fillnozzle (integer(IN) :: blockID,
!!                        real(IN)    :: dt,
!!                        real(IN)    :: time)
!!
!!  
!!  
!! DESCRIPTION
!!
!!  Apply the condition in the nozzle.
!!
!!
!! ARGUMENTS
!!
!!  blockID    : index of blocks to operate on
!!  dt         : current timestep
!!  time       : current time
!!
!!***

subroutine Heat_fillnozzle (blockID,dt,time)
!
!==============================================================================
!
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
    Grid_getCellCoords, Grid_putPointData, Grid_getDeltas
  !use Hydro_data, ONLY: hy_unsplitEosMode
  !use Eos_interface, ONLY : Eos_wrapped
  use Simulation_data
  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer,intent(IN) :: blockID
  real,intent(IN) :: dt,time

  integer :: i,j,k, istat
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  real :: dx, dy, dz
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  logical :: gcell = .true.
  real :: Br, Bz, Bphi
  real, pointer, dimension(:,:,:,:) :: solnData, &
        solnFaceXData, solnFaceYData, solnFaceZData

  integer :: nozzle=1
  real, dimension(3) :: cellvec, cellB, del
  real :: radius, length, sig, distance, theta, vel
  real, dimension(3) :: plnvec, jetvec, rvec, phivec

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(xCoord(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)
  call Grid_getDeltas(blockID, del)
  dx = del(1)
  dy = del(2)
  dz = del(3)
  
  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell, xCoord, sizeX)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell, zCoord, sizeZ)
  
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  call Grid_getBlkPtr(blockID,solnFaceXData,FACEX)
  call Grid_getBlkPtr(blockID,solnFaceYData,FACEY)
  call Grid_getBlkPtr(blockID,solnFaceZData,FACEZ)

  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH, IAXIS)
           cellvec = (/ xCoord(i), yCoord(j), zCoord(k) /)
           call hy_uhd_jetNozzleGeometry(nozzle,cellvec,radius,length,distance,sig,theta,jetvec,rvec,plnvec,phivec)
           if ( radius.lt.sim(nozzle)%radius.and.(abs(length).lt.sim(nozzle)%length) ) then
              ! inside the jet nozzle
              if (radius.gt.(sim(nozzle)%radius-sim(nozzle)%bfeather_outer)) then
                 vel = taper(nozzle, radius, sim(nozzle)%velocity, sim_windVel)
                 solnData(VELX_VAR,i,j,k) = vel*jetvec(1)*sig
                 solnData(VELY_VAR,i,j,k) = vel*jetvec(2)*sig
                 solnData(VELZ_VAR,i,j,k) = vel*jetvec(3)*sig
                 solnData(DENS_VAR,i,j,k) = taper(nozzle, radius, sim(nozzle)%density, sim_rhoAmbient)
                 solnData(PRES_VAR,i,j,k) = taper(nozzle, radius, sim(nozzle)%pressure, sim_pAmbient)
              else
                 solnData(VELX_VAR,i,j,k) = sim(nozzle)%velocity*jetvec(1)*sig
                 solnData(VELY_VAR,i,j,k) = sim(nozzle)%velocity*jetvec(2)*sig
                 solnData(VELZ_VAR,i,j,k) = sim(nozzle)%velocity*jetvec(3)*sig
                 solnData(DENS_VAR,i,j,k) = sim(nozzle)%density
                 solnData(PRES_VAR,i,j,k) = sim(nozzle)%pressure
              endif
              solnData(EINT_VAR,i,j,k) = solnData(PRES_VAR,i,j,k)/(solnData(GAME_VAR,i,j,k)-1.0)/solnData(DENS_VAR,i,j,k)
              !solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)+&
              !                      0.5*(solnData(VELX_VAR,i,j,k)**2 +&
              !                           solnData(VELY_VAR,i,j,k)**2 +&
              !                           solnData(VELZ_VAR,i,j,k)**2)
              solnData(JET_SPEC,i,j,k) = 1.0 - sim_smallX
              solnData(ISM_SPEC,i,j,k) = sim_smallX
              ! apply B field at cell center
              call hy_uhd_getBfield(nozzle,time,radius,length,0.0,Br,Bz,Bphi)
              cellB = Br*plnvec + Bz*jetvec + Bphi*phivec
              solnData(MAGX_VAR,i,j,k) = cellB(1)
              solnData(MAGY_VAR,i,j,k) = cellB(2)
              solnData(MAGZ_VAR,i,j,k) = cellB(3)
              solnData(MAGP_VAR,i,j,k) = sum(cellB*cellB)/(8.0*PI)
              ! update total energy
              solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)+&
                                    0.5*(solnData(VELX_VAR,i,j,k)**2 +&
                                         solnData(VELY_VAR,i,j,k)**2 +&
                                         solnData(VELZ_VAR,i,j,k)**2) +&
                                         sum(cellB*cellB)/(8.0*PI)
              
              ! apply B field at face center
              del = (/ -dx/2.0, 0.0, 0.0 /)
              call hy_uhd_jetNozzleGeometry(nozzle,cellvec+del,radius,length,distance,sig,theta,jetvec,rvec,plnvec,phivec)
              call hy_uhd_getBfield(nozzle,time,radius,length,0.0,Br,Bz,Bphi)
              cellB = Br*plnvec + Bz*jetvec + Bphi*phivec*sig
              solnFaceXData(MAG_FACE_VAR,i,j,k) = cellB(1)

              del = (/ 0.0, -dy/2.0, 0.0 /)
              call hy_uhd_jetNozzleGeometry(nozzle,cellvec+del,radius,length,distance,sig,theta,jetvec,rvec,plnvec,phivec)
              call hy_uhd_getBfield(nozzle,time,radius,length,0.0,Br,Bz,Bphi)
              cellB = Br*plnvec + Bz*jetvec + Bphi*phivec*sig
              solnFaceYData(MAG_FACE_VAR,i,j,k) = cellB(2)

              del = (/ 0.0, 0.0, -dz/2.0 /)
              call hy_uhd_jetNozzleGeometry(nozzle,cellvec+del,radius,length,distance,sig,theta,jetvec,rvec,plnvec,phivec)
              call hy_uhd_getBfield(nozzle,time,radius,length,0.0,Br,Bz,Bphi)
              cellB = Br*plnvec + Bz*jetvec + Bphi*phivec
              solnFaceZData(MAG_FACE_VAR,i,j,k) = cellB(3)

           endif
        enddo
     enddo
  enddo

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

  !call Eos_wrapped(hy_unsplitEosMode, blkLimits, blockID)

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  call Grid_releaseBlkPtr(blockID,solnFaceXData,FACEX)
  call Grid_releaseBlkPtr(blockID,solnFaceYData,FACEY)
  call Grid_releaseBlkPtr(blockID,solnFaceZData,FACEZ)

  return

end subroutine Heat_fillnozzle


