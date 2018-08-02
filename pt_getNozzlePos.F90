!!
!! DESCRIPTION
!!
!!  Get the random positions at the nozzle cross-section
!!  for jet particles.

subroutine pt_getNozzlePos(nAdd, pos)

  use Particles_data, ONLY : pt_meshMe
  use Simulation_data, ONLY : sim, sim_smallx, cross, sim_lowerRefHalf
  use Driver_data, ONLY : dr_dt

  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Flash_mpi.h"
#include "Particles.h"

  integer, INTENT(IN)   :: nAdd
  real, dimension(MDIM,nAdd), INTENT(OUT) :: pos
  real, dimension(MDIM) :: rxvec, ryvec
  ! The thikness of the layer for adding particles at the nozzle cross-section
  real          :: del
  real          :: r, theta, z, fac
  integer       :: i, nozzle=1, ierr

  rxvec = cross(sim(nozzle)%jetvec, (/ 0.0, sim_smallx, 1.0 /))
  rxvec = rxvec / sqrt(sum(rxvec(:)*rxvec(:)))
  ryvec = cross(sim(nozzle)%jetvec, rxvec)

  del = sim(nozzle)%velocity*dr_dt

  do i = 1, nAdd
     call RANDOM_NUMBER(r)
     r = sqrt(r)*sim(nozzle)%radius*0.8

     call RANDOM_NUMBER(theta)
     ! take square root of r to ensure uniform random distribution on the disk
     theta = theta*2.0*PI

     call RANDOM_NUMBER(z)
     ! z is a random number between 0 and 1
     ! We want lower half of the domain have (2**sim_lowerRefHalf)**3 times fewer
     ! particles to allow similar # of particles / block

     fac = (2**sim_lowerRefHalf)**MDIM
     z = z - 1.0/(fac + 1.0)

     if (z .lt. 0.0) then
         z = -sim(nozzle)%length + z*(fac+1.0)*del
     else
         z = sim(nozzle)%length + z*(fac+1.0)/fac*del
     endif

     pos(:,i) = sim(nozzle)%pos &
                + sim(nozzle)%jetvec*z&
                + r*(rxvec*cos(theta) + ryvec*sin(theta))
  enddo

end subroutine pt_getNozzlePos
