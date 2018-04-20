!! This subroutine reset the shock variables for shock.
!!
!! ARGUMENTS
!!
!!   particle : particle on which to operate
!!   shock_ind: index for shock
!!   dsa_ind  : diffusive shock acceleration power-law index

subroutine pt_resetShockVars(particle, shock_ind, dsa_ind)

   implicit none

#include "Flash.h"

   real,dimension(NPART_PROPS),intent(INOUT) :: particle
   integer, INTENT(in) :: shock_ind
   real, INTENT(in) :: dsa_ind

   particle(IND1_PART_PROP+shock_ind-1) = dsa_ind
   particle(TAU1_PART_PROP+shock_ind-1) = 1E-100
   particle(ICT1_PART_PROP+shock_ind-1) = 1E-100
   particle(DEN1_PART_PROP+shock_ind-1) = particle(DENS_PART_PROP)
   particle(AGE1_PART_PROP+shock_ind-1) = 0.0
   particle(WHCH_PART_PROP) = REAL(shock_ind) + 0.1

end subroutine pt_resetShockVars

