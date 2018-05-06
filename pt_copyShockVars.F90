!! This subroutine copies the shock variables from orig to dest.

subroutine pt_copyShockVars(particle, orig, dest)

   implicit none

#include "Flash.h"

   real,dimension(NPART_PROPS),intent(INOUT) :: particle
   integer, INTENT(in) :: dest, orig

   particle(IND1_PART_PROP+dest-1) = particle(IND1_PART_PROP+orig-1)
   particle(TAU1_PART_PROP+dest-1) = particle(TAU1_PART_PROP+orig-1)
   particle(DEN1_PART_PROP+dest-1) = particle(DEN1_PART_PROP+orig-1)
   particle(TAD1_PART_PROP+dest-1) = particle(TAD1_PART_PROP+orig-1)
   particle(ICT1_PART_PROP+dest-1) = particle(ICT1_PART_PROP+orig-1)

end subroutine pt_copyShockVars

