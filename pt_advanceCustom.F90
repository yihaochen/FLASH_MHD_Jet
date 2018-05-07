!!****if* source/Particles/localAPI/pt_advanceCustom
!!
!! NAME
!!
!!  pt_advanceCharged
!!
!! SYNOPSIS
!!
!!  call pt_advanceCharged(real(in)   :: dtOld,
!!                        real(in)    :: dtNew,
!!                        real(inout) :: particles(:,p_count),
!!                        integer(in) :: p_count)
!!                        integer(in) :: ind)
!!
!! DESCRIPTION
!!
!!   Advances particles in time, based on rungekutta (pt_advanceRK.F90).
!!
!! ARGUMENTS
!!
!!   dtOld : previous time interval 
!!   dtNew : current time interval
!!   particle -- particles on which to operate
!!   p_count - the number of particles in the list to advance
!!   ind -- index into pt_typeInfo
!!
!! NOTES
!!
!!  No special handling is done for the first call - it is assumed that particle
!!  initialization fills in initial velocity components properly.
!!
!!
!!***

subroutine pt_advanceCustom(dtOld,dtNew, particles,p_count, ind)
    
  
  use Particles_data, ONLY: pt_numLocal, pt_maxPerProc, &
       useParticles, pt_typeInfo, &
       pt_gcMaskForAdvance, pt_gcMaskSizeForAdvance, pt_meshMe, &
       pt_posAttrib, pt_velNumAttrib,pt_velAttrib, pt_meshNumProcs,&
       pt_customNumAttrib, pt_customAttrib

  use Driver_data, ONLY : dr_simTime, dr_initialSimTime, dr_globalMe
  use Hydro_data, ONLY : hy_bref
  use Grid_interface, ONLY : Grid_mapMeshToParticles
  use Simulation_data, ONLY : sim_ptAddPeriod, sim, sim_rCore
  
  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Flash_mpi.h"
#include "Particles.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  integer, INTENT(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(INOUT) :: particles

  integer       :: i,particleTypes
  
  real          :: jumpx,jumpy,jumpz
  real,allocatable :: origVel(:,:)
  integer :: part_props=NPART_PROPS
  integer :: mapType 
  
  integer      :: nozzle=1
  real         :: rho13, A, Aic, prob, dsa_ind, r2
!!------------------------------------------------------------------------------
  
  !write(*,'(i5, A28, i5)') pt_meshMe, '[pt_advance], pt_numLocal=', pt_numLocal
  !write(*,'(i5, A28, i5)') pt_meshMe, '[pt_advance], p_count    =', p_count
  ! Don't do anything if runtime parameter isn't set
  if (.not.useParticles .or. p_count.eq.0) return

  mapType=pt_typeInfo(PART_MAPMETHOD,ind)

  ! Particle initialization
  do i = 1, p_count
  ! if the particle is newly added, map the properties first
     if (particles(DEN0_PART_PROP,i) .le. 0.0 ) then
        ! This will update DEN0 and DENS
        call Grid_mapMeshToParticles(particles(:,i:i),&
             part_props,BLK_PART_PROP, 1,&
             pt_posAttrib,pt_customNumAttrib,pt_customAttrib,mapType)
        particles(DEN0_PART_PROP,i) = particles(DENS_PART_PROP,i)

        call Grid_mapMeshToParticles(particles(:,i:i),&
             part_props,BLK_PART_PROP, 1,&
             pt_posAttrib,pt_velNumAttrib,pt_velAttrib,QUADRATIC)

     endif
  enddo

  !!------------------------------------------------------------------------------
  ! Update the particle positions to temporary ("predicted") values
  do i = 1, p_count
 
     jumpx = dtNew * particles(VELX_PART_PROP,i)
     particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + jumpx

     if(NDIM >1) then
        jumpy = dtNew * particles(VELY_PART_PROP,i)
        particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + jumpy
        if(NDIM >2) then
           jumpz = dtNew * particles(VELZ_PART_PROP,i)
           particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + jumpz
        endif
     endif

  enddo


  ! Now save the original velocity values
  allocate(origVel(p_count,MDIM))
  origVel(:,1) = particles(VELX_PART_PROP,1:p_count)
  origVel(:,2) = particles(VELY_PART_PROP,1:p_count)
  origVel(:,3) = particles(VELZ_PART_PROP,1:p_count)

  ! Map the updated gas velocity field at the temporary positions to
  ! obtain a second estimate of velocities;

  call Grid_mapMeshToParticles(particles,&
       part_props, BLK_PART_PROP, p_count,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,QUADRATIC)

  ! Adjust particle positions, using the second point velocities
  do i = 1, p_count
     particles(POSX_PART_PROP,i) =  particles(POSX_PART_PROP,i) + &
          dtNew * 0.5*(particles(VELX_PART_PROP,i) - origVel(i,1))
     if(NDIM>1)&
          particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + &
          dtNew * 0.5*(particles(VELY_PART_PROP,i) - origVel(i,2))
     if(NDIM>2)&
          particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + &
          dtNew * 0.5*(particles(VELZ_PART_PROP,i) - origVel(i,3))
  enddo

  ! done with this temporary storage
  deallocate(origVel)



  ! Map the updated gas velocity field onto the current particle positions to
  ! obtain the updated particle velocities - for the next integration step
  ! as well as for particle plot files etc.

  call Grid_mapMeshToParticles(particles,&
       part_props, BLK_PART_PROP,p_count,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,QUADRATIC)

  ! Map the custom fields onto the partitcles
  call Grid_mapMeshToParticles(particles,&
       part_props,BLK_PART_PROP, p_count,&
       pt_posAttrib,pt_customNumAttrib,pt_customAttrib,mapType)

  ! update the synchrotron lifetime and cutoff gamma
  do i = 1, p_count
    if (particles(SHKS_PART_PROP,i) .gt. 1.0) then

       ! Energy power-law index for diffusive shock acceleration
       ! SHKS is the compression ratio
       dsa_ind = (particles(SHKS_PART_PROP,i)+2.0) / (particles(SHKS_PART_PROP,i)-1.0)
       if (dsa_ind .lt. particles(IND1_PART_PROP,i)) then
          ! This particle is in a shock front that is stronger than it
          ! encountered before, we reset the cooling integral

          ! If it is not already in SHOCK1, copy existing shocks first
          if (particles(WHCH_PART_PROP,i) .gt. 3.0) then
             call pt_copyShockVars(particles(:,i), 2, 3)
             call pt_copyShockVars(particles(:,i), 1, 2)
          else if (particles(WHCH_PART_PROP,i) .gt. 2.0) then
             ! This will copy the current cooling history to IND2
             call pt_copyShockVars(particles(:,i), 1, 2)
          endif
          call pt_resetShockVars(particles(:,i), 1, dsa_ind, dr_simTime)
          particles(GAMC_PART_PROP,i) = 1E100
       else if ( (dsa_ind .lt. particles(IND2_PART_PROP,i)) .and.&
                 (particles(WHCH_PART_PROP,i) .gt. 2.0) ) then
          ! Shock is stronger than previously stored IND2 and not in the
          ! remaining end of IND1

          ! If it is not already in IND2, copy IND2 to IND3 first
          if (particles(WHCH_PART_PROP,i) .gt. 3.0) then
             ! This will copy the current cooling history to IND2
             call pt_copyShockVars(particles(:,i), 2, 3)
          endif
          ! This particle is in a shok, but the shock strength is weaker than
          ! the strongest shock it encountered
          call pt_resetShockVars(particles(:,i), 2, dsa_ind, dr_simTime)
       else if ( (dsa_ind .lt. particles(IND3_PART_PROP,i)) .and.&
                 (particles(WHCH_PART_PROP,i) .gt. 3.0) ) then
          call pt_resetShockVars(particles(:,i), 3, dsa_ind, dr_simTime)
       endif
    else
       ! Outside of a shock
       particles(WHCH_PART_PROP,i) = 100.0
    endif
    ! A = 4/3 * sigmaT * c * beta^2 / (me*c^2) * U
    ! sigmaT = 8*pi / 3 * e^4 / (c^4*me^2)
    ! A = 32*pi / 9 * e^4/(me^3*c^5) * U
    !   = 4/9 * e^4/(me^3*c^5) * Bcgs^2
    !
    ! e^4/(me^3*c^5) = 2.907728E-9
    ! if units == "none" (default) -> hy_bref=1.0
    ! if units == "cgs" -> hy_bref=sqrt(4*pi)
    A = 32.0*PI/9.0*2.907728E-9/hy_bref/hy_bref/2.0*&
        sum(particles(MAGX_PART_PROP:MAGZ_PART_PROP,i)*particles(MAGX_PART_PROP:MAGZ_PART_PROP,i))
    ! U_cmb = 4.17E-13 erg/cm^3 *(1+z)^4  using T = 2.725 K
    Aic = 32.0*PI/9.0*2.907728E-9*4.17E-13
    rho13 = (particles(DENS_PART_PROP,i)/particles(DEN0_PART_PROP,i))**(1.0/3.0)
    ! Distance square to the jet origin, assuming to be the center of the galaxy
    r2 = sum( (particles(POSX_PART_PROP:POSZ_PART_PROP,i)-sim(nozzle)%pos)**2 )
    ! Cooling integration for the injection tracer
    if (particles(DEN0_PART_PROP,i) .gt. 0.0) then
       ! Synchrotron cooling
       particles(TAU0_PART_PROP,i) = particles(TAU0_PART_PROP,i) + rho13*A*dtNew
       ! Inverse-Compton cooling of stellar light
       particles(ICT0_PART_PROP,i) = particles(ICT0_PART_PROP,i) + rho13*Aic/(1+r2/sim_rCore**2)*dtNew
       ! Inverse-Compton cooling of CMB
       particles(CMBT_PART_PROP,i) = particles(CMBT_PART_PROP,i) + rho13*Aic*dtNew
    endif
    ! Cooling integration for the 1st shock tracer
    if ( (particles(DEN1_PART_PROP,i) .gt. 0.0) .and.&
         (abs(particles(WHCH_PART_PROP,i)-1.0) .gt. 0.1) ) then
       particles(TAU1_PART_PROP,i) = particles(TAU1_PART_PROP,i) + rho13*A*dtNew
       particles(ICT1_PART_PROP,i) = particles(ICT1_PART_PROP,i) + rho13*Aic/(1+r2/sim_rCore**2)*dtNew
       particles(GAMC_PART_PROP,i) = rho13 / particles(TAU1_PART_PROP,i)
    endif
    ! Cooling integration for the 2nd shock tracer
    if ( (particles(DEN2_PART_PROP,i) .gt. 0.0) .and.&
         (abs(particles(WHCH_PART_PROP,i)-2.0) .gt. 0.1) ) then
       particles(TAU2_PART_PROP,i) = particles(TAU2_PART_PROP,i) + rho13*A*dtNew
       particles(ICT2_PART_PROP,i) = particles(ICT2_PART_PROP,i) + rho13*Aic/(1+r2/sim_rCore**2)*dtNew
    endif
    ! Cooling integration for the 3rd shock tracer
    if ( (particles(DEN3_PART_PROP,i) .gt. 0.0) .and.&
         (abs(particles(WHCH_PART_PROP,i)-3.0) .gt. 0.1) ) then
       particles(TAU3_PART_PROP,i) = particles(TAU3_PART_PROP,i) + rho13*A*dtNew
       particles(ICT3_PART_PROP,i) = particles(ICT3_PART_PROP,i) + rho13*Aic/(1+r2/sim_rCore**2)*dtNew
    endif

  enddo

    !write(*,'(i5, A28, i5)') pt_meshMe, '[pt_advance_end], pt_numLocal=', pt_numLocal
    !write(*,'(i5, A28, i5)') pt_meshMe, '[pt_advance_end], p_count    =', p_count
  
end subroutine pt_advanceCustom
