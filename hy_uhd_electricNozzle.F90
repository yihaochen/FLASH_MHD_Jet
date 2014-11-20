Subroutine hy_uhd_electricNozzle(blockID, blkLimits, blkLimitsGC)

  use Simulation_data
  use Driver_data,      ONLY : dr_simTime,dr_dt,dr_globalMe
  use Grid_interface,   ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr,Grid_getCellCoords,&
                               Grid_getDeltas

  implicit none
  
#include "constants.h"
#include "Flash.h"

  integer, intent(IN)  :: blockID
  integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
  integer :: i,j,k,xyz, istat
  integer :: sizeX,sizeY,sizeZ
  logical :: gcell = .true.
  integer :: nozzle = 1
  real, pointer, dimension(:,:,:,:) :: E

  real, dimension(3) :: nozvec, nozvecf, edgevec, advect, torvec, polvec
  real, dimension(3) :: jetvec, rvec, plnvec, phivec, phivecOld
  real :: Ar, Az, Aphi, Arold, Azold, Aphiold!, thetavel, angvel, vel

  real :: length,radius,distance, theta, cellsig!, fac, fillfac
  real :: Efac!, torfac!, blockingRTaper, toroidalReplenishLTaper

  ! For debug
  real, dimension(3) :: del
  real :: dx, dy, dz
  call Grid_getDeltas(blockID, del)
  dx = del(1)
  dy = del(2)
  dz = del(3)
  
  
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(sim_xCoord(sizeX),stat=istat)
  allocate(sim_yCoord(sizeY),stat=istat)
  allocate(sim_zCoord(sizeZ),stat=istat)
  allocate(sim_xCoordf(sizeX+1),stat=istat)
  allocate(sim_yCoordf(sizeY+1),stat=istat)
  allocate(sim_zCoordf(sizeZ+1),stat=istat)

  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell, sim_xCoord, sizeX)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell, sim_yCoord, sizeY)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell, sim_zCoord, sizeZ)
  call Grid_getCellCoords(IAXIS,blockID,FACES,gcell, sim_xcoordf, sizeX+1)
  call Grid_getCellCoords(JAXIS,blockID,FACES,gcell, sim_ycoordf, sizeX+1)
  call Grid_getCellCoords(KAXIS,blockID,FACES,gcell, sim_zcoordf, sizeX+1)

  call Grid_getBlkPtr(blockID,E,SCRATCH)

  !write(*,*) blockID
  
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
   do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
    do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
        !write (*,*) 'blockID, ijk, xyz=', blockID,i,j,k,xyz
        ! Vector to cell center
        nozvec(:)=[sim_xcoord(i),sim_ycoord(j),sim_zcoord(k)]
        !write(*,*) 'nozvec =', nozvec
        
        ! Vector to cell corner (where faces meet)
        nozvecf(:)=[sim_xcoordf(i),sim_ycoordf(j),sim_zcoordf(k)]
        !write(*,*) 'nozvecf =', nozvecf


     ! Now loop over the three faces
     do xyz=EX_SCRATCH_GRID_VAR, EZ_SCRATCH_GRID_VAR
        ! edge-centered vector
        edgevec(:)   = nozvecf(:)
        edgevec(xyz) = nozvec(xyz)
        !write(*,*) 'edgevec =', edgevec

        call hy_uhd_jetNozzleGeometry(nozzle,edgevec,radius,length, &
             distance,cellsig,theta,jetvec,rvec,plnvec,phivec)

        ! now calculate the additional vector potential to be added
        
        !call hy_uhd_jetNozzleFill(rvec,jetvec,radius,length,distance,cellsig,&
        !     vel,dens,pres,eint,fac,fillfac )
        call hy_uhd_getA(nozzle, radius, length, 0.0, Ar, Az, Aphi)
        
        !
        ! Taper factors - smoothly transition from imposed to FLASH solution
        !
        !LRTaper = 0.5*(1.0+cos(PI*max(0.0,(min(1.0,&
        !          (abs(length)-sim(nozzle)%length)/sim(nozzle)%zfeather)))))*&
        !          taperR(nozzle, radius, 1.0, 0.0)
        Efac = ETaper(nozzle, radius, length, 0.0, 0.0, 1.0)
        
        !torfac = taperL(nozzle, length, 0.0, 1.0)
        
        ! nozzle face taper factor for toroidal field
        ! 0 means use injection scheme, 
        ! 1 means use FLASH calculated value
        
        !blockingLTaper = 0.5*(1.0+cos(PI*max(-1.0,(min(1.0,&
        !                 (abs(length)-sim(nozzle)%zTorInj)/sim(nozzle)%zfeather)))))
        
        
        ! radial taper factor
        !blockingRTaper = taperR(nozzle, radius, 1.0, 0.0)

        ! replacement for missing v x B. This is the z-derivative 
        ! of blockingLTaper
        !toroidalReplenishLTaper=0.5*PI / sim(nozzle)%zfeather * &
        !     sin(PI * (max(-1.0,min(1.0 , &
        !     (abs(length) - sim(nozzle)%zTorInj) / sim(nozzle)%zfeather ))))

        ! rotation velicity to be used for advection
        !thetavel(:)=cross(sim(nozzle)%angvel(:),rvec(:)*distance)
        
        
        ! inject or advect?
        if ( dr_simTime.ge.sim(nozzle)%timeMHDon .and. &
            (dr_simTime-dr_dt).lt.sim(nozzle)%timeMHDon ) then
           
           ! First jet timestep - inject initial field

           !if (dr_globalMe==MASTER_PE) then
           !   write(*,*) '        Initial EM field injection'
           !endif
           advect(:)=[0.,0.,0.]
           
           ! add field
           E(xyz,i,j,k) = E(xyz,i,j,k) - (Ar*plnvec(xyz)+Az*jetvec(xyz))&
                                       - Aphi*phivec(xyz)/dr_dt

        else if ((dr_simTime-dr_dt).ge.sim(nozzle)%timeMHDon) then
           
           ! Toroidal field update: Replace the flux advected away

           ! Az has dimension [B*length/time]
           ! We don't need to divide it by dr_dt again.
           torvec(xyz) = - (Ar*plnvec(xyz)+Az*jetvec(xyz))
           
           !call hy_uhd_jetNozzleGeometryOld(edgevec,radius,length, &
           !     distance,cellsig,theta,jetvec,rvec,plnvec,phivec)
           !call hy_uhd_getA(radius,length,distance,phi,theta,&
           !     jetvec,rvec,plnvec,phivec,aphiold,azold,aold)
           

           call hy_uhd_jetNozzleGeometryOld(nozzle,edgevec,radius,length, &
                distance,cellsig,theta,jetvec,rvec,plnvec,phivecOld)
           call hy_uhd_getAOld(nozzle, radius, length, 0.0, Arold, Azold, Aphiold)
           ! rotate and move poloidal field (subtract old field, add moved)
           advect(xyz) = -(Aphi*phivec(xyz) - Aphiold*phivecOld(xyz))/dr_dt

           ! Finally, update the E-field so that FLASH can update B using
           ! dB = - (curl E) * dt
           ! ( dA/dt = -E )
           
           ! The mixing factors are included in torvec and advect.
           E(xyz,i,j,k) = E(xyz,i,j,k)*Efac + torvec(xyz) + advect(xyz)
           
        endif

     enddo ! End of the loop over 3 directions
    enddo
   enddo
  enddo

  deallocate(sim_xCoord)
  deallocate(sim_yCoord)
  deallocate(sim_zCoord)
  deallocate(sim_xCoordf)
  deallocate(sim_yCoordf)
  deallocate(sim_zCoordf)
  call Grid_releaseBlkPtr(blockID,E,SCRATCH)

contains

  function ETaper(nozzle, r, z, var_cen, var_in, var_out)
  ! Taper function for electric field in both R and z direction
  ! The electric nozzle is bigger than the hydro nozzle to avoid boundary
  ! magnification of the magnetic field due to the large velocity.
  ! See taper() in Simulation_data.F90 for more details.

    integer, INTENT(in) :: nozzle
    real, INTENT(in) :: r, z, var_cen, var_in, var_out
    real :: r1, r2, rout, zout, zjet, var_zcen, var_zin
    real :: ETaper
    r1 = sim(nozzle)%rFeatherIn
    r2 = sim(nozzle)%radius + sim(nozzle)%rFeatherOut
    rout = r2 + sim(nozzle)%rFeatherOut
    zjet = sim(nozzle)%length
    zout = zjet + sim(nozzle)%zFeather

    ! z part
    if (abs(z).ge.0.0 .and. abs(z).lt.zjet) then
      var_zin = 1.0*var_in
      var_zcen = 1.0*var_cen
    else if (abs(z).ge.zjet .and. abs(z).lt.zout) then
      var_zin = (-2*(zout-abs(z))**3/(zout-zjet)**3 &
      + 3*(zout-abs(z))**2/(zout-zjet)**2)*(var_in-var_out) + var_out
      var_zcen = (-2*(zout-abs(z))**3/(zout-zjet)**3 &
      + 3*(zout-abs(z))**2/(zout-zjet)**2)*(var_cen-var_out) + var_out
    else
      var_zin = var_out
      var_zcen = var_out
    end if

    ! radial part
    if (r.ge.0.0 .and. r.lt.r1) then
      ETaper = (-2*r**3/r1**3 + 3*r**2/r1**2)*(var_zin-var_zcen) + var_zcen
    else if (r.ge.r1 .and. r.lt.r2) then
      ETaper = 1.0*var_zin
    else if (r.ge.r2 .and. r.lt.rout) then
      ETaper = (-2*(rout-r)**3/(rout-r2)**3 &
      + 3*(rout-r)**2/(rout-r2)**2)*(var_zin-var_out) + var_out
    else
      ETaper = var_out
    endif

  end function ETaper

End Subroutine hy_uhd_electricNozzle
