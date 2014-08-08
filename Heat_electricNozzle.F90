Subroutine Heat_electricNozzle( nozzle,i,j,k,E )

  use Simulation_data
  use Driver_data,      ONLY : dr_simTime,dr_nstep,dr_dt
  use constants

  implicit none
  
#include "constants.h"
#include "Flash.h"

  integer, intent(IN) :: i,j,k, nozzle
  real, dimension(3), intent(INOUT) :: E
  integer :: xyz

  real, dimension(3) :: nozvec, nozvecf, edgevec, advect, torvec, polvec
  real, dimension(3) :: jetvec, rvec, plnvec, phivec
  real :: Ar, Az, Aphi, Arold, Azold, Aphiold!, thetavel, angvel, vel

  real :: length,radius,distance, theta, cellsig!, fac, fillfac
  real :: LRTaper!blockingLTaper, blockingRTaper, toroidalReplenishLTaper
  
  ! now calculate the additional vector potential to be added
  
  !
  !   Nozzle geometry:
  !
  !               ______________
  !              /              \
  !             /                \
  !             |                |
  !             |\______________/|
  !             |  _____________ |
  !             | /             \|
  !             |/               |
  !             |                |
  !              \______________/
  !
  
  ! Vector to cell center
  nozvec(:)=[sim_xcoord(i),sim_ycoord(j),sim_zcoord(k)]
  
  ! Vector to cell corner (where faces meet)
  nozvecf(:)=[sim_xcoordf(i),sim_ycoordf(j),sim_zcoordf(k)]
  
  ! Now loop over the three faces
  do xyz=EX_SCRATCH_GRID_VAR, EZ_SCRATCH_GRID_VAR

     ! edge-centered vector
     edgevec(:)=nozvecf(:)
     edgevec(xyz)=nozvec(xyz)
     
     !write(*,*) '###########################################################1'
     call hy_uhd_jetNozzleGeometry(nozzle,edgevec,radius,length, &
          distance,cellsig,theta,jetvec,rvec,plnvec,phivec)
     !write(*,*) '###########################################################2'
     
     
     !call hy_uhd_jetNozzleFill(rvec,jetvec,radius,length,distance,cellsig,&
     !     vel,dens,pres,eint,fac,fillfac )
     call hy_uhd_getA(nozzle, dr_simTime, radius, length, 0.0, Ar, Az, Aphi)
     call hy_uhd_getA(nozzle, dr_simTime-dr_dt, radius, length, 0.0, Arold, Azold, Aphiold)
     
     !
     ! Taper factors - smoothly transition from imposed to FLASH solution
     !
     LRTaper = 0.5*(1.0+cos(PI*max(-1.0,(min(1.0,&
               (abs(length)-sim(nozzle)%bPosZ)/sim(nozzle)%zfeather)))))
     
     ! nozzle face taper factor for toroidal field
     ! 1 means use injection scheme, 
     ! 0 means use FLASH calculated value
     
     !blockingLTaper = 0.5*(1.0+cos(PI*max(-1.0,(min(1.0,&
     !                 (abs(length)-sim(nozzle)%bPosZ)/sim(nozzle)%zfeather)))))
     
     
     ! radial taper factor
     !blockingRTaper = taperR(nozzle, radius, 1.0, 0.0)

     ! replacement for missing v x B. This is the z-derivative 
     ! of blockingLTaper
     !toroidalReplenishLTaper=0.5*PI / sim(nozzle)%zfeather * &
     !     sin(PI * (max(-1.0,min(1.0 , &
     !     (abs(length) - sim(nozzle)%bPosZ) / sim(nozzle)%zfeather ))))

     ! rotation velicity to be used for advection
     !thetavel(:)=cross(sim(nozzle)%angvel(:),rvec(:)*distance)
     
     
     ! inject or advect?
     if (dr_simTime.ge.sim(nozzle)%timeMHDon .and. &
          dr_simTime.lt.sim(nozzle)%timeMHDon+dr_dt) then
        
        ! First jet timestep - inject initial field
        
        ! Vector potential for poloidal field is 
        ! phihat*poloidalGeneratingFunction*distance*sin(theta)/2
        
        advect(:)=[0.,0.,0.]
        
        ! add field
        !if (i .eq. 8 .and. j.eq.8 .and. k.eq.8) then
        !write(*,*) 'Injet initial field'
        !endif
        E(xyz) = E(xyz) - (Ar*plnvec(xyz)+Az*jetvec(xyz))&!/dr_dt&
                 !*sim(nozzle)%zfeather/sim(nozzle)%velocity &
                 -Aphi*phivec(xyz)/dr_dt

     else if (dr_simTime.ge.sim(nozzle)%timeMHDon+dr_dt) then
        
        ! Toroidal field update: Replace the flux advected away
        !torvec(:)= - Ar(:)*sum(jetvec(:)*sim(nozzle)%velocity)*&
        !     toroidalReplenishLTaper
        torvec(xyz) = - (Ar*plnvec(xyz)+Az*jetvec(xyz))!*blockingLTaper
        
        !call hy_uhd_jetNozzleGeometryOld(edgevec,radius,length, &
        !     distance,cellsig,theta,jetvec,rvec,plnvec,phivec)
        !call hy_uhd_getA(radius,length,distance,phi,theta,&
        !     jetvec,rvec,plnvec,phivec,aphiold,azold,aold)
        
        ! rotate and move poloidal field (subtract old field, add moved)
        advect(xyz) = -(Aphi - Aphiold)*phivec(xyz)/dr_dt

        ! Finally, update the E-field so that FLASH can update B using
        ! dB = - (curl E) * dt
        
        !if (i .eq. 8 .and. j.eq.8 .and. k.eq.8) then
        !write(*,*) 'Advect field'
        !endif
        E(xyz) = &
             E(xyz)*(1.0 - LRTaper)+&!1.0*blockingRTaper) + &
             advect(xyz) + torvec(xyz)
        
        if (i .eq. 8 .and. j.eq.8 .and. k.eq.8) then
           !write(*,*) 'Ar', Ar, torvec(:)
           !write(*,*)'E',xyz,E(xyz)
        endif
        
     endif

  enddo

End Subroutine Heat_electricNozzle
