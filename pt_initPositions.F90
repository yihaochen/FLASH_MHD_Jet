!!****if* source/Simulation/SimulationMain/magnetoHD/MHD_Jet/
!!
!! NAME
!!    pt_initPositions
!!
!! SYNOPSIS
!!
!!    call pt_initPositions(integer(in)  :: blockID,
!!                          logical(OUT) :: success)
!!
!! DESCRIPTION
!!    Initialize particle locations.  This version sets up particles
!!      which are evenly distributed in circular space
!!
!! ARGUMENTS
!!
!!  blockID:        local block ID containing particles to create
!!  success:        returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!! PARAMETERS
!!
!!    pt_numX:      number of particles along physical x-axis of domain
!!    pt_numY:      number of particles along physical y-axis of domain
!!    pt_numZ:      number of particles along physical z-axis of domain
!!
!!***


subroutine pt_initPositions (blockID,success)

  implicit none

  integer, INTENT(in) :: blockID
  logical,intent(OUT) :: success

  success=.true.

  return

!----------------------------------------------------------------------

end subroutine pt_initPositions


