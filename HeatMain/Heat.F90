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
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getDeltas
  use hy_uhd_interface, ONLY : hy_uhd_staggeredDivb
  !use Hydro_data, ONLY: hy_unsplitEosMode
  !use Eos_interface, ONLY : Eos_wrapped
  use Driver_data, ONLY : dr_globalMe, dr_nStep
  use Simulation_data
  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN)::blockList
  real,intent(IN) :: dt,time
  real, dimension(MDIM) :: del
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  logical :: halfTimeAdvance = .false.

  integer :: blockID, blkInd, nozzle=1

  !call calc_jet(nozzle, time)
  !if (dr_globalMe==MASTER_PE .and. mod(dr_nStep,20)==0) then
  !   write(*,'(a,2es11.3, f7.2)') '      (p, rho, M)=', &
  !   sim(nozzle)%pressure, sim(nozzle)%density, sim(nozzle)%mach
  !endif
  
  do blkInd=1,blockCount
     blockID = blockList(blkInd)
     call Grid_getDeltas(blockID,del)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     call Heat_fillnozzle(blockID,dt,time)

     !call hy_uhd_staggeredDivb(blockID,dt,del,blkLimits,blkLimitsGC,halfTimeAdvance)

  enddo

  return

end subroutine Heat


