!!***** source/Simulation/SimulationMain/magnetoHD/MHD_Jet/Simulation_jiggle.F90
!!
!! NAME
!!
!!  Simulation_jiggle
!!
!! SYNOPSIS
!!
!!  Simulation_jiggle(  integer (IN) :: nozzle,
!!                      real    (IN) :: time,
!!                      real    (IN) :: dt      )
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!  nozzle       -  index of the nozzle (currently 1)
!!  time         -  current simulation time
!!  dt           -  current simulation time step
!!
!!***

subroutine Simulation_jiggle( nozzle, time, dt )

    use Driver_data, ONLY: dr_globalMe
    use Simulation_data

    implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

    !! ---- Argument List ----------------------------------
    integer, INTENT(IN) :: nozzle
    real, INTENT(in) :: time, dt
    !! -----------------------------------------------------

    real, dimension(3), save :: ex, ey, ez, exprime, eyprime, ezprime, v, vprime, r, rcart
    real, dimension(3) :: phihat,thetahat, jetaxis
    real, save :: thetajet, phijet
    real :: vphi, vtheta, dv, sigma, rn, psi, dummyn
    integer :: ierr
    logical, save :: first_call=.true.
    real :: svr, radius, f, wigglerad
    real, dimension(3) :: rvec, vr

    !! End of data declaration ***********************************************

    if (first_call) then
       !
       ! Note:
       !
       ! sim(nozzle)%jetvec is the jet vector in simulation coordinates
       !
       ! thetajet and phijet are the polar angles of the jet axis in
       ! precession cone coordinates
       !

       ! cartesian vectors
       ex = [1.0, 0.0, 0.0]
       ey = [0.0, 1.0, 0.0]
       ez = [0.0, 0.0, 1.0]

       ! base vectors for precession cone
       ezprime(:)=sim(nozzle)%coneVec(:)

       if (ezprime(1).eq.ez(1).and.ezprime(2).eq.ez(2).and.ezprime(3).eq.ez(3)) then
          exprime(:)=ex(:)
       else
          exprime(:)=cross(ezprime(:),-ey(:))
          exprime(:)=exprime/sqrt(sum(exprime(:)**2))
       endif

       eyprime(:)=cross(ezprime(:), exprime(:))
       eyprime(:)=eyprime/sqrt(sum(eyprime(:)**2))

       rcart(:)=sim(nozzle)%jetvec(:)

       ! jet axis in precession cone coordinates
       r(:)=[sum(rcart(:)*exprime(:)),&
            sum(rcart(:)*eyprime(:)),&
            sum(rcart(:)*ezprime(:))]
       r(:)=r(:)/sqrt(sum(r(:)**2))

       ! jet angles in precession cone base
       thetajet=acos(r(3))
       if (abs(r(2)).gt.0.or.abs(r(1)).gt.0) then
          phijet=2.0*atan(r(2)/(r(1) + &
               sqrt(r(1)**2 + r(2)**2)))
       else
          phijet=PI
       endif

       ! base vectors in precession cone coordinates
       phihat(:)=cross(ez(:),r(:))
       if (sum(phihat(:)**2).gt.0) then
          phihat(:)=phihat(:)/sqrt(sum(phihat(:)**2))
       else
          phihat(:)=ey(:)
       endif

       thetahat(:)=cross(phihat(:),r(:))

       vprime(:) = sim(nozzle)%jetvec(:) - sim(nozzle)%jetvecOld(:)
       ! velocity in precession cone coordinates
       v(:) = [sum(exprime*vprime), sum(eyprime*vprime), sum(ezprime*vprime)]
       ! v is always unity
       if (sum(v**2).gt.0) then
          v(:)=v(:)/sqrt(sum(v(:)**2))
       endif
       !if (dr_globalMe == MASTER_PE) write (*,*) 'v = ', v
       !if (dr_globalMe == MASTER_PE) write (*,*) 'agV=', sim(nozzle)%angVel
       !if (dr_globalMe == MASTER_PE) write (*,*) 'v = ', v

       !sim(nozzle)%jetvec(:)=ezprime(:)
       !sim(nozzle)%jetvec_old(:)=sim(nozzle)%jetvec(:)
       first_call=.false.

    endif


    ! initializing?
    !if (time.ge.sim(nozzle)%tOn .and. &
    !     time.lt.(sim(nozzle)%tOn + dt)) then
    !
    !   ! Initialize jet vector
    !   sim(nozzle)%switch=.true.

    if (time.ge.sim(nozzle)%tOn .and. &
        time.lt.sim(nozzle)%tOn+sim(nozzle)%duration) then

       !if (dr_globalMe==MASTER_PE) then
       !   write(*,'(A28, 3f8.4)') '[Simulation_jiggle] jetvec  ', sim(nozzle)%jetvec
       !   write(*,'(A28, 3f8.4)') '[Simulation_jiggle] rcart   ', rcart
       !   write(*,'(A28, 3f8.4)') '[Simulation_jiggle] r       ', r
       !   write(*,'(A28, 2f8.4)') '[Simulation_jiggle] th, phi ', thetajet, phijet
       !   write(*,'(A28, 3f8.4)') '[Simulation_jiggle] thetahat', thetahat
       !   write(*,'(A28, 3f8.4)') '[Simulation_jiggle] phihat  ', phihat
       !endif

    !if (MyPE.eq.MasterPE) then
    !   write(*,*)'jiggle  1: ',sim(nozzle)%theta,sim(nozzle)%phi
    !   write(*,*)'jiggle  2: ',sim(nozzle)%pos
    !endif


       ! Initialize jet vector
       !sim(nozzle)%switch=.true.

       phihat(:)=cross(ez(:),r(:))
       if (sqrt(sum(phihat(:)**2)).gt.0) then
          phihat=phihat/sqrt(sum(phihat(:)**2))
       else
          phihat(:)=ey(:)
       endif

       thetahat(:)=cross(phihat(:),r(:))

       ! Now express velocity in new coordinate system
       vphi=sum(v(:)*phihat(:))
       vtheta=sum(v(:)*thetahat(:))

       ! Random walk step
       sigma=2.0*PI*exp(-(thetajet/sim(nozzle)%precangle)**4)
       call RANDOM_SEED(put=sim(nozzle)%randSeed)
       call RANDOM_NUMBER(rn)
       call MPI_Bcast(rn,1,MPI_DOUBLE_PRECISION,MASTER_PE,MPI_COMM_WORLD,ierr)
       call RANDOM_SEED(get=sim(nozzle)%randSeed)
       ! taper precession cone
       psi=(rn - 0.5)*sigma - PI/2.0

       dv = sim(nozzle)%nutation*sqrt(dt/sim(nozzle)%duration)
       ! Add to velocity
       vphi=vphi + dv*cos(psi)
       vtheta=vtheta + dv*sin(psi)

       !dummyn=sqrt(vphi**2 + vtheta**2)
       !if (dummyn.gt.0.0) then
       !   vphi=vphi/dummyn
       !   vtheta=vtheta/dummyn
       !endif

       ! New velocity in coordinates relative to cone axis
       v(:)=vphi*phihat(:) + vtheta*thetahat(:)

       ! And normalize
       v(:)=v(:)/sqrt(sum(v(:)**2))

       !if (dr_globalMe == MASTER_PE) write (*,*) 'v = ', v
       ! Now move jet axis and normalize
       r(:)=r(:) + v(:)*sim(nozzle)%precession*dt/sim(nozzle)%duration
       r(:)=r(:)/sqrt(sum(r(:)**2))
       !if (dr_globalMe == MASTER_PE) write (*,*) 'r = ', r

       ! Finally, calculate new angles
       thetajet=acos(r(3))
       !if (abs(r(1)).gt.0.or.abs(r(2)).gt.0) then
       !   phijet=2.0*atan(r(2)/&
       !        (r(1) + sqrt(r(1)**2 + r(2)**2)))
       !else
       !   phijet=PI
       !endif

       ! Now express in simulation coordinates
       !jprime(:)=[cos(phijet)*sin(thetajet),&
       !     sin(phijet)*sin(thetajet),cos(thetajet)]
       !if (dr_globalMe == MASTER_PE) write (*,*) 'jp= ', jprime
       !jetaxis(:)=exprime(:)*jprime(1) + &
       !     eyprime(:)*jprime(2) + ezprime(:)*jprime(3)
       jetaxis(:)=exprime(:)*r(1) + eyprime(:)*r(2) + ezprime(:)*r(3)

       ! sim(nozzle)%theta and sim(nozzle)%phi are measured in simulations coordinates
       !sim(nozzle)%theta=acos(jetaxis(3)/sqrt(sum(jetaxis(:)**2)))
       !if (jetaxis(2).gt.-sqrt(jetaxis(1)**2 + jetaxis(2)**2)) then
       !   sim(nozzle)%phi=2.0*atan(jetaxis(2)/(jetaxis(1) + &
       !        sqrt(jetaxis(1)**2 + jetaxis(2)**2)))
       !else
       !   sim(nozzle)%phi=pi
       !endif

       ! angular velocity for precession
       sim(nozzle)%angVel(:)=&
            cross(sim(nozzle)%jetvec(:),jetaxis(:))/dt
       !if (dr_globalMe == MASTER_PE) write (*,*) 'agV=', sim(nozzle)%angVel

       ! velocity in simulation coordinates (for consistancey at restart)
       vprime(:) = jetaxis(:)-sim(nozzle)%jetvec(:)
       ! velocity in precession cone coordinates
       v(:) = [sum(exprime*vprime), sum(eyprime*vprime), sum(ezprime*vprime)]
       v(:)=v(:)/sqrt(sum(v(:)**2))
       !if (dr_globalMe == MASTER_PE) write (*,*) 'v" = ', v

       ! move this to jiggle
       sim(nozzle)%jetvecOld(:)=sim(nozzle)%jetvec(:)
       sim(nozzle)%jetvec(:)=jetaxis(:)

       !if (dr_globalMe==MASTER_PE) then
       !   write(*,'(A28, 2f9.5)') '[Simulation_jiggle] angle : ',thetajet, phijet
       !   !write(*,'(A28, 3f9.5)') '[Simulation_jiggle] vector: ',sim(nozzle)%jetvec
       !endif

    !else
    !   sim(nozzle)%switch=.false.
    endif


end subroutine Simulation_jiggle

