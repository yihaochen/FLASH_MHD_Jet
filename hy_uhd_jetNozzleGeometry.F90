Subroutine hy_uhd_jetNozzleGeometry(nozzle,cellvec,radius,length, &
     distance,sig,theta,jetvec,rvec,plnvec,phivec )

  use Simulation_data

  implicit none
  
  integer, intent(IN) :: nozzle
  real, dimension(3), intent(IN) :: cellvec
  real, intent(OUT) :: radius, length, sig, distance, theta
  real, dimension(3), intent(OUT) :: plnvec, jetvec, rvec, phivec
  real, dimension(3) :: nozvec
  real :: norm

  ! Jer direction
  jetvec(:)=sim(nozzle)%jetvec

  ! Vector to cell center
  nozvec(:)=cellvec(:) - sim(nozzle)%pos(:)

  norm=sqrt(sum(nozvec(:)**2))
  if (norm.gt.0) then
     rvec=nozvec/norm
  else
     rvec=jetvec
  endif

  ! Position along the jet
  length=sum(jetvec(:)*nozvec(:))

  ! Left or right of nozzle?
  sig=1.0
  if (length.lt.0) sig=-1.0

  ! Polar angle
  distance=sqrt(sum(nozvec(:)*nozvec(:)))
  if (distance.gt.0) then
     theta=acos(length/distance)
  else
     theta=0.
  endif
  
  ! Vector from jet axis to cell center
  plnvec=nozvec(:)-jetvec(:)*length
 
  ! Radial distance from jet axis
  radius=sqrt(sum(plnvec(:)*plnvec(:)))
  
  ! Normalized radius vector
  if (radius.ne.0.) then 
     plnvec(:)=plnvec(:)/radius
  else
     plnvec(:)=[0.,0.,0.]
  endif

  ! phivector
  phivec(:)=cross(jetvec,rvec)
  norm=sum(phivec(:)*phivec(:))
  if (norm.gt.0) then
     phivec(:)=phivec(:)/sqrt(norm)
  else
     phivec(:)=[jetvec(2),-jetvec(1),0]
  endif

End Subroutine hy_uhd_jetNozzleGeometry
