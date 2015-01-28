!!****f* source/Simulation/Simulation_sendOutputData
!!
!! NAME
!!  Simulation_sendOutputData
!!
!! SYNOPSIS
!! 
!!  Simulation_sendOutputData()
!!  
!! DESCRIPTION 
!!
!! This routine sends the scalar variables owned by the Simulation unit
!! to the IO unit, to be written to a checkpoint file.
!!
!!
!!***

subroutine Simulation_sendOutputData()

  use Simulation_data
  use IO_interface, ONLY :  IO_setScalar

  implicit none
  
#include "Simulation.h"
  
  integer :: nozzle = 1

  call IO_setScalar('nozzlePosX', sim(nozzle)%pos(1))
  call IO_setScalar('nozzlePosY', sim(nozzle)%pos(2))
  call IO_setScalar('nozzlePosZ', sim(nozzle)%pos(3))
  call IO_setScalar('nozzleVecX', sim(nozzle)%jetvec(1))
  call IO_setScalar('nozzleVecY', sim(nozzle)%jetvec(2))
  call IO_setScalar('nozzleVecZ', sim(nozzle)%jetvec(3))
  call IO_setScalar('nozzleAngVelX', sim(nozzle)%angVel(1))
  call IO_setScalar('nozzleAngVelY', sim(nozzle)%angVel(2))
  call IO_setScalar('nozzleAngVelZ', sim(nozzle)%angVel(3))
  call IO_setScalar('nozzleLinVelX', sim(nozzle)%linVel(1))
  call IO_setScalar('nozzleLinVelY', sim(nozzle)%linVel(2))
  call IO_setScalar('nozzleLinVelZ', sim(nozzle)%linVel(3))
  !call IO_setScalar('nozzlet0', sim(nozzle)%t0)
  

end subroutine Simulation_sendOutputData

