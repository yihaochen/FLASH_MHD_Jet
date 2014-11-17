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
  !real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  real :: dx, dy, dz
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  logical :: gcell = .true.
  real :: Br, Bz, Bphi
  real, pointer, dimension(:,:,:,:) :: solnData

  integer :: nozzle=1
  real, dimension(3) :: cellvec, del
  real :: radius, length, sig, distance, theta, vel, fac
  real, dimension(3) :: plnvec, jetvec, rvec, phivec, velvec

  if (present(init_in)) then
      init = init_in
  endif

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(sim_xCoord(sizeX),stat=istat)
  allocate(sim_yCoord(sizeY),stat=istat)
  allocate(sim_zCoord(sizeZ),stat=istat)
  !allocate(sim_xCoordf(sizeX+1),stat=istat)
  !allocate(sim_yCoordf(sizeY+1),stat=istat)
  !allocate(sim_zCoordf(sizeZ+1),stat=istat)
  call Grid_getDeltas(blockID, del)
  dx = del(1)
  dy = del(2)
  dz = del(3)
  
  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell, sim_xCoord, sizeX)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell, sim_yCoord, sizeY)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell, sim_zCoord, sizeZ)

  !call Grid_getCellCoords(IAXIS,blockID,FACES,gcell, sim_xcoordf, sizeX+1)
  !call Grid_getCellCoords(JAXIS,blockID,FACES,gcell, sim_ycoordf, sizeX+1)
  !call Grid_getCellCoords(KAXIS,blockID,FACES,gcell, sim_zcoordf, sizeX+1)
  
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  !call Grid_getBlkPtr(blockID,solnFaceXData,FACEX)
  !call Grid_getBlkPtr(blockID,solnFaceYData,FACEY)
  !call Grid_getBlkPtr(blockID,solnFaceZData,FACEZ)
  !call Grid_getBlkPtr(blockID,E,SCRATCH)

  ! Set the electric field to 0 to avoid duplicated advection when 
  ! hy_uhd_staggeredDivB is called again in Heat.F90
  !E(EX_SCRATCH_GRID_VAR:EZ_SCRATCH_GRID_VAR,:,:,:) = 0.0


  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
   do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
    do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
       cellvec = (/ sim_xCoord(i), sim_yCoord(j), sim_zCoord(k) /)
       call hy_uhd_jetNozzleGeometry(nozzle,cellvec,radius,length,distance,&
                                     sig,theta,jetvec,rvec,plnvec,phivec)

       if ((radius.le.(sim(nozzle)%radius+sim(nozzle)%rFeatherOut))&
           .and.(abs(length).le.(sim(nozzle)%length+sim(nozzle)%zFeather))) then
       ! inside the jet nozzle
          vel = sim(nozzle)%velocity*sin(PI/2.0*min(abs(length),sim(nozzle)%length)*sig/sim(nozzle)%length)
          if ((radius.le.sim(nozzle)%radius) .and. (abs(length).le.sim(nozzle)%length)) then
              fac = 1.0
          else
              fac = 0.0
          endif
          !fac = taper(nozzle, radius, length, 1.0, 1.0, 0.0)
          ! smooth transition from nozzle to flash grid
          ! 1.0 for nozzle injection; 0.0 for flash solution
          vel = sim(nozzle)%velocity*&
                sin(PI/2.0*min(abs(length),sim(nozzle)%length)*sig/sim(nozzle)%length)
          velvec = vel*jetvec&
                   + sim(nozzle)%outflowR*sim(nozzle)%velocity*plnvec*&
                     0.5*(1.0+cos(PI*(min(0.0, radius-sim(nozzle)%radius)/sim(nozzle)%rFeatherOut)))&
                   + sim(nozzle)%linVel*fac + cross(sim(nozzle)%angVel,rvec*distance)
          solnData(VELX_VAR:VELZ_VAR,i,j,k) = velvec*fac + solnData(VELX_VAR:VELZ_VAR,i,j,k)*(1.0-fac)
          solnData(DENS_VAR,i,j,k) = sim(nozzle)%density*fac+&
                                     solnData(DENS_VAR,i,j,k)*(1.0-fac)
          solnData(PRES_VAR,i,j,k) = sim(nozzle)%pressure*fac+&
                                     solnData(PRES_VAR,i,j,k)*(1.0-fac)

          solnData(EINT_VAR,i,j,k) = solnData(PRES_VAR,i,j,k)/solnData(DENS_VAR,i,j,k)&
                                     /(solnData(GAME_VAR,i,j,k)-1.0)
          solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)+&
                                0.5*(solnData(VELX_VAR,i,j,k)**2 +&
                                     solnData(VELY_VAR,i,j,k)**2 +&
                                     solnData(VELZ_VAR,i,j,k)**2)
          solnData(JET_SPEC,i,j,k) = (fac - sim_smallX) + solnData(JET_SPEC,i,j,k)*(1.0-fac)
          solnData(ISM_SPEC,i,j,k) = sim_smallX +  solnData(ISM_SPEC,i,j,k)*(1.0-fac)
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

          
          !if (i .eq. 8 .and. j.eq.8 .and. k.eq.8) then
          !write(*,*)'EX',E(EX_SCRATCH_GRID_VAR,i,j,k),fac
          !write(*,*)'EY',E(EY_SCRATCH_GRID_VAR,i,j,k)
          !write(*,*)'EZ',E(EZ_SCRATCH_GRID_VAR,i,j,k)
          !endif

       endif ! inside the nozzle

       ! Fill the nozzle with toroidal and poloidal field and close the field
       ! outside the nozzle
       !call Heat_electricNozzle(nozzle,i,j,k,E(:,i,j,k))
    enddo
   enddo
  enddo

  deallocate(sim_xCoord)
  deallocate(sim_yCoord)
  deallocate(sim_zCoord)
  !deallocate(sim_xCoordf)
  !deallocate(sim_yCoordf)
  !deallocate(sim_zCoordf)

  !call Eos_wrapped(hy_unsplitEosMode, blkLimits, blockID)

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  !call Grid_releaseBlkPtr(blockID,solnFaceXData,FACEX)
  !call Grid_releaseBlkPtr(blockID,solnFaceYData,FACEY)
  !call Grid_releaseBlkPtr(blockID,solnFaceZData,FACEZ)
  !call Grid_releaseBlkPtr(blockID,E,SCRATCH)

  return

end subroutine Heat_fillnozzle


