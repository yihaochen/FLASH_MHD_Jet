

subroutine pt_getRandomPos(pos)

  use Particles_data, ONLY : pt_meshMe
  use Simulation_data, ONLY : sim, sim_smallx, cross

  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Particles.h"

  !real, dimension(MDIM), INTENT(IN) :: del
  real :: del=0.01
  real, dimension(MDIM), INTENT(OUT) :: pos
  real, dimension(MDIM) :: rxvec, ryvec
  real          :: r, theta, z
  integer       :: nozzle=1

  rxvec = cross(sim(nozzle)%jetvec, (/ 0.0, sim_smallx, 1.0 /))
  rxvec = rxvec / sqrt(sum(rxvec(:)*rxvec(:)))
  ryvec = cross(sim(nozzle)%jetvec, rxvec)

  call RANDOM_NUMBER(r)
  call RANDOM_NUMBER(theta)
  call RANDOM_NUMBER(z)
  
  ! take square root of r to ensure uniform random distribution on the disk
  r = sqrt(r)*(sim(nozzle)%radius+sim(nozzle)%rFeatherOut)
  theta = theta*2.0*PI
  z = z - 0.5
  if (z .lt. 0.0) then
      z = (-1.0+z*del)*sim(nozzle)%length
  else
      z = ( 1.0+z*del)*sim(nozzle)%length
  endif

  pos = sim(nozzle)%pos &
        + sim(nozzle)%jetvec*z&
        + r*(rxvec*cos(theta) + ryvec*sin(theta))

end subroutine pt_getRandomPos
