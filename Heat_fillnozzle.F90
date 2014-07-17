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

subroutine Heat_fillnozzle (blockID,dt,time,init_in)
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
  logical,intent(IN),optional :: init_in
  logical :: init = .false.

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
  real :: radius, length, sig, distance, theta, vel, fac
  real, dimension(3) :: plnvec, jetvec, rvec, phivec

  if (present(init_in)) then
      init = init_in
  endif

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
  !call Grid_getBlkPtr(blockID,solnFaceXData,FACEX)
  !call Grid_getBlkPtr(blockID,solnFaceYData,FACEY)
  !call Grid_getBlkPtr(blockID,solnFaceZData,FACEZ)

  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
   do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
       cellvec = (/ xCoord(i), yCoord(j), zCoord(k) /)
       call hy_uhd_jetNozzleGeometry(nozzle,cellvec,radius,length,distance,sig,theta,jetvec,rvec,plnvec,phivec)
       if (init) then
          !if ((radius.le.2.0*(sim(nozzle)%radius+sim(nozzle)%bfeather_outer))&
          !    .and.(abs(length).le.2.0*(sim(nozzle)%length+sim(nozzle)%zfeather))) then
          !   if ((radius.gt.2.0*sim(nozzle)%radius).or.(abs(length).gt.2.0*sim(nozzle)%length)) then
          !      solnData(DENS_VAR,i,j,k) = taper(nozzle, 0.5*radius, 0.5*length, sim(nozzle)%density,&
          !                                       sim(nozzle)%density, sim_rhoAmbient)
          !      solnData(PRES_VAR,i,j,k) = taper(nozzle, 0.5*radius, 0.5*length, sim(nozzle)%pressure,&
          !                                       sim(nozzle)%pressure, sim_pAmbient)
          !   else
          !      solnData(DENS_VAR,i,j,k) = sim(nozzle)%density
          !      solnData(PRES_VAR,i,j,k) = sim(nozzle)%pressure
          !   endif
          !   solnData(EINT_VAR,i,j,k) = solnData(PRES_VAR,i,j,k)/solnData(DENS_VAR,i,j,k)&
          !                              /(solnData(GAME_VAR,i,j,k)-1.0)
          !   solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)+&
          !                         0.5*(solnData(VELX_VAR,i,j,k)**2 +&
          !                              solnData(VELY_VAR,i,j,k)**2 +&
          !                              solnData(VELZ_VAR,i,j,k)**2)
          !   solnData(JET_SPEC,i,j,k) = 1.0 - sim_smallX
          !   solnData(ISM_SPEC,i,j,k) = sim_smallX
          !endif
            

          !if (init) then
          !! inside the jet nozzle radius, but outside the jet nozzle length
          !! only for smooth initial condition
          !    if ((abs(length).le.(sim(nozzle)%length+sim(nozzle)%zfeather)).and. &
          !        (abs(length).gt.sim(nozzle)%length)) then
          !       vel = taper_zout(nozzle, radius, length, sim(nozzle)%velocity, sim(nozzle)%velocity, sim_windVel)
          !       solnData(VELX_VAR,i,j,k) = vel*jetvec(1)*sig
          !       solnData(VELY_VAR,i,j,k) = vel*jetvec(2)*sig
          !       solnData(VELZ_VAR,i,j,k) = vel*jetvec(3)*sig
          !       solnData(DENS_VAR,i,j,k) = taper_zout(nozzle, radius, length,&
          !                                  sim(nozzle)%density, sim(nozzle)%density, sim_rhoAmbient)
          !       solnData(PRES_VAR,i,j,k) = taper_zout(nozzle, radius, length,&
          !                                  sim(nozzle)%pressure, sim(nozzle)%pressure, sim_pAmbient)
          !       solnData(EINT_VAR,i,j,k) = solnData(PRES_VAR,i,j,k)/(solnData(GAME_VAR,i,j,k)-1.0)/solnData(DENS_VAR,i,j,k)
          !       solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)+&
          !                             0.5*(solnData(VELX_VAR,i,j,k)**2 +&
          !                                  solnData(VELY_VAR,i,j,k)**2 +&
          !                                  solnData(VELZ_VAR,i,j,k)**2)
          !       call hy_uhd_getBfield(nozzle,time,radius,length,0.0,Br,Bz,Bphi)
          !       cellB = Br*plnvec + Bz*jetvec + Bphi*phivec
          !       solnData(MAGX_VAR,i,j,k) = cellB(1)
          !       solnData(MAGY_VAR,i,j,k) = cellB(2)
          !       solnData(MAGZ_VAR,i,j,k) = cellB(3)
          !       solnData(MAGP_VAR,i,j,k) = sum(cellB*cellB)/(8.0*PI)
          !       solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)+&
          !                             0.5*(solnData(VELX_VAR,i,j,k)**2 +&
          !                                  solnData(VELY_VAR,i,j,k)**2 +&
          !                                  solnData(VELZ_VAR,i,j,k)**2)+&
          !                                  sum(cellB*cellB)/(8.0*PI)
          !       ! apply B field at face center
          !       del = (/ -dx/2.0, 0.0, 0.0 /)
          !       call hy_uhd_jetNozzleGeometry(nozzle,cellvec+del,radius,length,distance,sig,theta,jetvec,rvec,plnvec,phivec)
          !       call hy_uhd_getBfield(nozzle,time,radius,length,0.0,Br,Bz,Bphi)
          !       cellB = Br*plnvec + Bz*jetvec + Bphi*phivec
          !       solnFaceXData(MAG_FACE_VAR,i,j,k) = cellB(1)

          !       del = (/ 0.0, -dy/2.0, 0.0 /)
          !       call hy_uhd_jetNozzleGeometry(nozzle,cellvec+del,radius,length,distance,sig,theta,jetvec,rvec,plnvec,phivec)
          !       call hy_uhd_getBfield(nozzle,time,radius,length,0.0,Br,Bz,Bphi)
          !       cellB = Br*plnvec + Bz*jetvec + Bphi*phivec
          !       solnFaceYData(MAG_FACE_VAR,i,j,k) = cellB(2)

          !       del = (/ 0.0, 0.0, -dz/2.0 /)
          !       call hy_uhd_jetNozzleGeometry(nozzle,cellvec+del,radius,length,distance,sig,theta,jetvec,rvec,plnvec,phivec)
          !       call hy_uhd_getBfield(nozzle,time,radius,length,0.0,Br,Bz,Bphi)
          !       cellB = Br*plnvec + Bz*jetvec + Bphi*phivec
          !       solnFaceZData(MAG_FACE_VAR,i,j,k) = cellB(3)
          !    endif
          !endif

       endif

       if ((radius.le.sim(nozzle)%radius+sim(nozzle)%bfeather_outer)&
           .and.(abs(length).le.sim(nozzle)%length+sim(nozzle)%zfeather)) then
       ! inside the jet nozzle
          fac = taperR(nozzle, radius, 0.0, 1.0)
          ! smooth transition from nozzle to flash grid
          ! 0.0 for nozzle injection; 1.0 for flash solution
          solnData(VELX_VAR,i,j,k) = sim(nozzle)%velocity*jetvec(1)*sig*(1.0-fac)+&
                                     solnData(VELX_VAR,i,j,k)*fac
          solnData(VELY_VAR,i,j,k) = sim(nozzle)%velocity*jetvec(2)*sig*(1.0-fac)+&
                                     solnData(VELY_VAR,i,j,k)*fac
          solnData(VELZ_VAR,i,j,k) = sim(nozzle)%velocity*jetvec(3)*sig*(1.0-fac)+&
                                     solnData(VELZ_VAR,i,j,k)*fac
          solnData(DENS_VAR,i,j,k) = sim(nozzle)%density*(1.0-fac)+&
                                     solnData(DENS_VAR,i,j,k)*fac
          solnData(PRES_VAR,i,j,k) = sim(nozzle)%pressure*(1.0-fac)+&
                                     solnData(PRES_VAR,i,j,k)*fac

          solnData(EINT_VAR,i,j,k) = solnData(PRES_VAR,i,j,k)/solnData(DENS_VAR,i,j,k)&
                                     /(solnData(GAME_VAR,i,j,k)-1.0)
          solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)+&
                                0.5*(solnData(VELX_VAR,i,j,k)**2 +&
                                     solnData(VELY_VAR,i,j,k)**2 +&
                                     solnData(VELZ_VAR,i,j,k)**2)
          solnData(JET_SPEC,i,j,k) = (1.0 - sim_smallX - fac) + solnData(JET_SPEC,i,j,k)*fac
          solnData(ISM_SPEC,i,j,k) = sim_smallX + solnData(ISM_SPEC,i,j,k)*fac
          !! apply B field at cell center
          !call hy_uhd_getBfield(nozzle,time,radius,length,0.0,Br,Bz,Bphi)
          !cellB = Br*plnvec + Bz*jetvec + Bphi*phivec
          !solnData(MAGX_VAR,i,j,k) = cellB(1)
          !solnData(MAGY_VAR,i,j,k) = cellB(2)
          !solnData(MAGZ_VAR,i,j,k) = cellB(3)
          !solnData(MAGP_VAR,i,j,k) = sum(cellB*cellB)/(8.0*PI)
          !! update total energy
          !solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)+&
          !                      0.5*(solnData(VELX_VAR,i,j,k)**2 +&
          !                           solnData(VELY_VAR,i,j,k)**2 +&
          !                           solnData(VELZ_VAR,i,j,k)**2) +&
          !                           sum(cellB*cellB)/(8.0*PI)
          !
          !! apply B field at face center
          !del = (/ -dx/2.0, 0.0, 0.0 /)
          !call hy_uhd_jetNozzleGeometry(nozzle,cellvec+del,radius,length,distance,sig,theta,jetvec,rvec,plnvec,phivec)
          !call hy_uhd_getBfield(nozzle,time,radius,length,0.0,Br,Bz,Bphi)
          !cellB = Br*plnvec + Bz*jetvec + Bphi*phivec
          !solnFaceXData(MAG_FACE_VAR,i,j,k) = cellB(1)

          !del = (/ 0.0, -dy/2.0, 0.0 /)
          !call hy_uhd_jetNozzleGeometry(nozzle,cellvec+del,radius,length,distance,sig,theta,jetvec,rvec,plnvec,phivec)
          !call hy_uhd_getBfield(nozzle,time,radius,length,0.0,Br,Bz,Bphi)
          !cellB = Br*plnvec + Bz*jetvec + Bphi*phivec
          !solnFaceYData(MAG_FACE_VAR,i,j,k) = cellB(2)

          !del = (/ 0.0, 0.0, -dz/2.0 /)
          !call hy_uhd_jetNozzleGeometry(nozzle,cellvec+del,radius,length,distance,sig,theta,jetvec,rvec,plnvec,phivec)
          !call hy_uhd_getBfield(nozzle,time,radius,length,0.0,Br,Bz,Bphi)
          !cellB = Br*plnvec + Bz*jetvec + Bphi*phivec
          !solnFaceZData(MAG_FACE_VAR,i,j,k) = cellB(3)

       endif ! inside the nozzle
    enddo
   enddo
  enddo

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

  !call Eos_wrapped(hy_unsplitEosMode, blkLimits, blockID)

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  !call Grid_releaseBlkPtr(blockID,solnFaceXData,FACEX)
  !call Grid_releaseBlkPtr(blockID,solnFaceYData,FACEY)
  !call Grid_releaseBlkPtr(blockID,solnFaceZData,FACEZ)

  return

end subroutine Heat_fillnozzle


