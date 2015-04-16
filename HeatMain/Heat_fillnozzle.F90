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
  use Eos_interface, ONLY : Eos_wrapped
  use Simulation_data
  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer,intent(IN) :: blockID
  real,intent(IN) :: dt,time

  integer :: i,j,k, istat
  integer,dimension(3) :: minijk
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  logical :: gcell = .true.
  real :: Br, Bz, Bphi
  real, pointer, dimension(:,:,:,:) :: solnData
  real, dimension(NPROP_VARS) :: mixData

  integer :: nozzle=1
  real, dimension(3) :: cellvec, del, mixvec
  real :: radius, length, sig, distance, theta, vel, facR, facL
  real, dimension(3) :: plnvec, jetvec, rvec, phivec, voutvec, velvec
  real :: r2, bf, rout, rmix

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(sim_xCoord(sizeX),stat=istat)
  allocate(sim_yCoord(sizeY),stat=istat)
  allocate(sim_zCoord(sizeZ),stat=istat)
  call Grid_getDeltas(blockID, del)
  
  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell, sim_xCoord, sizeX)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell, sim_yCoord, sizeY)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell, sim_zCoord, sizeZ)

  call Grid_getBlkPtr(blockID,solnData,CENTER)

  bf = sim(nozzle)%rFeatherOut
  r2 = sim(nozzle)%radius
  rout = r2 + bf
  rmix = rout + sim(nozzle)%rFeatherMix;

  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
   do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
       cellvec = (/ sim_xCoord(i), sim_yCoord(j), sim_zCoord(k) /)
       call hy_uhd_jetNozzleGeometry(nozzle,cellvec,radius,length,distance,&
                                     sig,theta,jetvec,rvec,plnvec,phivec)
       if ((radius.le.rmix).and.&
           (abs(length).le.(sim(nozzle)%length+sim(nozzle)%zFeatherMix))) then
          !! inside the jet nozzle
          !!if ((radius.le.rout) .and. (abs(length).le.sim(nozzle)%length)) then
          !!    fac = 1.0
          !!else
          !!    fac = 0.0
          !!endif
          facR = taperR(nozzle, radius, 1.0, 0.0, zero_center=.false.)
          facL = taperL(nozzle, length, 1.0, 0.0)
          ! smooth transition from nozzle to flash grid
          ! 1.0 for nozzle injection; 0.0 for flash solution

          ! vector to the outter boundary of the nozzle feather
          mixvec = cellvec+plnvec*(rmix-radius)
          ! mixData contains the values of the variables for interpolation
          call ht_getValueAtPoint(blockID, mixvec, del, mixData)

          vel = sim(nozzle)%velocity&
                *( 0.5*(1.0+cos(PI*(max(0.0, min(1.0, (radius-r2)/bf)))))&
                *(1.0-sim(nozzle)%outflowR)+sim(nozzle)%outflowR ) &
                *sin(PI/2.0*min(abs(length),sim(nozzle)%length)*sig/sim(nozzle)%length)
          voutvec = sim(nozzle)%outflowR*sim(nozzle)%velocity*plnvec&
                    !*coshat(radius-0.5*(r2+2.0*bf), 0.5*(r2+bf), bf, 1.0)
                    *0.5*(1.0+cos(PI*( max(-1.0, min(0.0,(radius-rout)/bf)) )))

          velvec = vel*jetvec + voutvec &
                   + sim(nozzle)%linVel + cross(sim(nozzle)%angVel,rvec*distance)
          solnData(VELX_VAR:VELZ_VAR,i,j,k) = &
                   (velvec*facR + mixData(VELX_VAR:VELZ_VAR)*(1.0-facR))*facL&
                   + solnData(VELX_VAR:VELZ_VAR,i,j,k)*(1.0-facL)
          solnData(DENS_VAR,i,j,k) = max(sim_smlrho, &
                   (sim(nozzle)%density*facR + mixData(DENS_VAR)*(1.0-facR))*facL&
                   + solnData(DENS_VAR,i,j,k)*(1.0-facL) ) 
          solnData(PRES_VAR,i,j,k) = (sim(nozzle)%pressure*facR + mixData(PRES_VAR)*(1.0-facR))*facL&
                   + solnData(PRES_VAR,i,j,k)*(1.0-facL)

          solnData(EINT_VAR,i,j,k) = solnData(PRES_VAR,i,j,k)/solnData(DENS_VAR,i,j,k)&
                                     /(solnData(GAME_VAR,i,j,k)-1.0)
          solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)&
                                     + 0.5*(solnData(VELX_VAR,i,j,k)**2 &
                                          + solnData(VELY_VAR,i,j,k)**2 &
                                          + solnData(VELZ_VAR,i,j,k)**2)
          solnData(JET_SPEC,i,j,k) = max( sim_smallX, facR + mixData(JET_SPEC)*(1.0-facR) )
          solnData(ISM_SPEC,i,j,k) = max( sim_smallX, mixData(ISM_SPEC)*(1.0-facR) )

       endif ! inside the nozzle

       ! Fill the nozzle with toroidal and poloidal field and close the field
       ! outside the nozzle
       !call Heat_electricNozzle(nozzle,i,j,k,E(:,i,j,k))
    enddo
   enddo
  enddo
  if (minval(solnData(DENS_VAR,:,:,:)) < 0.0) then
     minijk = minloc(solnData(DENS_VAR,:,:,:))
     write(*,'(i4, a, e11.3, 3i4)') blockID, ' min dens =', minval(solnData(DENS_VAR,:,:,:)), minijk
     write(*,'(a, 3es11.3)') '(x,y,z) =', sim_xCoord(minijk(1)), sim_yCoord(minijk(2)), sim_zCoord(minijk(3))
  endif

  deallocate(sim_xCoord)
  deallocate(sim_yCoord)
  deallocate(sim_zCoord)

  call Eos_wrapped(MODE_DENS_PRES, blkLimits, blockID)
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  return

end subroutine Heat_fillnozzle


