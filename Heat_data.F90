!!****if* source/Heat/HeatMain/MHD_Jet/Heat_data
!!
!! NAME
!!  Heat_data
!!
!! SYNOPSIS
!!
!!  use Heat_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for the wind tunnel problem with a step
!!
!! ARGUMENTS
!!
!!
!!
!!
!!
!!***

module Heat_data

  implicit none

#include "constants.h"
#include "Simulation.h"

  integer,save :: nPtProc
  real,allocatable,save,dimension(:,:) ::  pos

end module Heat_data


