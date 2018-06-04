!!
!! DESCRIPTION
!!
!!  Get the random positions at the nozzle cross-section
!!  for jet particles.

subroutine pt_getRandomPos(nAdd, pos)

  use Particles_data, ONLY : pt_meshMe
  use Simulation_data, ONLY : sim, sim_smallx, cross, sim_onlyHalf

  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Flash_mpi.h"
#include "Particles.h"

  !real, dimension(MDIM), INTENT(IN) :: del
  ! Ratio of the thikness of the layer to the nozzle length
  ! at the nozzle cross-section
  real :: del=0.01
  integer, INTENT(IN)   :: nAdd
  real, dimension(nAdd,MDIM), INTENT(OUT) :: pos
  real, dimension(MDIM) :: rxvec, ryvec
  real          :: r, theta, z
  integer       :: i, nozzle=1, ierr

  rxvec = cross(sim(nozzle)%jetvec, (/ 0.0, sim_smallx, 1.0 /))
  rxvec = rxvec / sqrt(sum(rxvec(:)*rxvec(:)))
  ryvec = cross(sim(nozzle)%jetvec, rxvec)

  do i = 1, nAdd
     call RANDOM_NUMBER(r)
     call RANDOM_NUMBER(theta)
     call RANDOM_NUMBER(z)
     
     ! take square root of r to ensure uniform random distribution on the disk
     r = sqrt(r)*sim(nozzle)%radius*0.8
     theta = theta*2.0*PI
     if (sim_onlyHalf) then
        z = (1.0+0.5*z*del)*sim(nozzle)%length
     else
        z = z - 0.5
        if (z .lt. 0.0) then
            z = (-1.0+z*del)*sim(nozzle)%length
        else
            z = ( 1.0+z*del)*sim(nozzle)%length
        endif
     endif

     pos(i,:) = sim(nozzle)%pos &
                + sim(nozzle)%jetvec*z&
                + r*(rxvec*cos(theta) + ryvec*sin(theta))
  enddo

end subroutine pt_getRandomPos
