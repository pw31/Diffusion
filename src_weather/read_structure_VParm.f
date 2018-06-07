**********************************************************************
      SUBROUTINE READ_STRUCTURE_VPARM(longitude, latitude)
**********************************************************************
      use NATURE,ONLY: bar,bk,amu,km,pi,grav
      use PARAMETERS,ONLY: struc_file,Teff,logg,Hp,vzconst,pconst,beta
      use READMODEL,ONLY: Rlay,Tlay,play,rholay,glay,vconvlay,
     >                    zlay,mulay,Difflay,Nlayers
      implicit none
      integer,intent(IN) :: longitude, latitude
      integer,parameter :: maxLayers=1000
      integer :: i,j,it,CN,numlayer,long,lat,ibase
      integer :: up,low,NNN
      real*8 :: planet_mass,planet_radius,meanv,meanv2
      real*8 :: dum,rr,p,T,nH,rho,g,mu,wmix,err,pnorm,Tnorm
      real*8 :: R0,P0,T0,rho0,drho,z,dz,qbase,integral
      real*8 :: Hplay,grad,U,V,bp,bT,Wpa,Wm
      real*8 :: ngas,lmean,Kn,vth,Dmicro
      real*8,dimension(maxLayers) :: Ulay,Vlay,Wpalay,Wmlay,vz,vz0
      character(len=30) :: longi, lati
      logical,parameter :: newMethod=.true.

      planet_mass   = 1.979614D+31  ! [g]
      planet_radius = 8.328818D+9   ! [cm]

******************************
* Finding Coordinates Number *
******************************
**********************************************************************
* Coordinates number : 
* (-180,0) -> 1     (-180,45) -> 2
* (-135,0) -> 3     (-135,45) -> 4
* (-90,0) -> 5      (-90,45) -> 6
* (-45,0) -> 7      (-45,45) -> 8
* (0,0) -> 9        (0,45) -> 10
* (45,0) -> 11      (45,45) -> 12
* (90,0) -> 13      (90,45) -> 14
* (135,0) -> 15     (135,45) -> 16
**********************************************************************

      if (longitude.eq.-180) then
        CN = 1
      else if (longitude.eq.-135) then
        CN = 3
      else if (longitude.eq.-90) then
        CN = 5
      else if (longitude.eq.-45) then
        CN = 7
      else if (longitude.eq.0) then
        CN = 9
      else if (longitude.eq.45) then
        CN = 11
      else if (longitude.eq.90) then
        CN = 13
      else if (longitude.eq.135) then
        CN = 15
      endif
      if (latitude.eq.45) then
        CN = CN+1
      endif

      write(*,*) 
      write(*,*) "Longitude = ",longitude
      write(*,*) "Latitude = ",latitude
      write(*,*) "reading 2D+1 PTv structure "//trim(struc_file)
      write(*,*) "=============================================="
      open(42, file=struc_file,status='old')
      read(42,*)
      read(42,*) Teff
      write(*,*) "Teff = ",Teff
      read(42,*)
      read(42,*)
      read(42,*) Nlayers
      write(*,*) "Nlayers = ",Nlayers
      read(42,*)
      read(42,*)
      read(42,*)
      write(*,*) "Searching the column ..."
      do j=1,CN-1
         do i=1,Nlayers
            read(42,*) numlayer, long, lat, bp, bT, U, V, Wpa, Wm
         enddo
         read(42,*)
      enddo
      write(*,*) "Reading the column ..."
      do i=1,Nlayers
         read(42,*) numlayer, long, lat, play(i), Tlay(i), Ulay(i),
     &              Vlay(i), Wpalay(i), Wmlay(i)
      end do
      write(*,*) "Longitude =",long
      write(*,*) "Latitude =",lat
      close(42)
      write(*,*) "2D+1 (P,T,v) structure closed"
      do i=1,Nlayers
         Wmlay(i) = Wmlay(i)*100.d0
      end do
      do i=1,Nlayers
         play(i) = play(i)*bar
      end do

*     -------------------------
*     ***  calculate mulay  ***
*     -------------------------
* a constant value is used, assuming a H2-dominated gas
      write(*,*) "Calculating Mulay ..."
      do i=1,Nlayers
        mulay(i)=2.3
c       write(*,*) "mulay = ",mulay(i)
      end do
      write(*,*) "mulay calculated"

*     --------------------------
*     ***  calculate rholay  ***
*     --------------------------
* the Ideal Gas Law is used to calculate Rho

      write(*,*) "Calculating Rholay ..."
      do i=1,Nlayers
        rholay(i) = mulay(i)*amu*play(i)/(bk*Tlay(i))
c       write(*,*) "rholay(",i,") = ",rholay(i)
      end do
      write(*,*) "Rholay calculated"

*     ----------------------------------
*     ***  calculate glay and log g  ***
*     ----------------------------------
* Constant value used, using the radius and the mass of the planet

      write(*,*) "Calculating glay ..."
      do i=1,Nlayers
        glay(i) = grav*planet_mass/planet_radius**2
c       write(*,*) "glay = ",glay(i)
      end do

      logg = DLOG10(glay(1))
      write(*,*) "glay and logg calculated"

*     ------------------------
*     ***  calculate Rlay  ***
*     ------------------------
* To calculate the height of each layer, the hydrostatic equilibrium is 
* used.  A reference point R0=radius of the planet is taken, with 
* P0=1bar, and it iterates from this point to each limits

      write(*,*) "Calculating Rlay ..."
      R0 = planet_radius  !Planet Radius, cm
      P0 = play(38)       !dyn cm-2
      T0 = Tlay(38)       !layer 38 is the closest one to 1 bar, so I take 
                          !the temperature of this layer as T0
      rho0 = mulay(38)*amu*P0/(bk*T0)     !rho calculated with the ideal gas law
      Rlay(38) = R0
      z = 0.0
      do i=37,1,-1 !the iteration begins from the reference point, 
                   !and go to the upper limit
         drho = sqrt(rholay(i)*rholay(i+1))
         z = z + ((play(i+1)-play(i))/(glay(i)*drho))
         Rlay(i) = R0+z
      end do
      z = 0.0
      do i=39,Nlayers              !another iteration goes to the lower limit
         drho = sqrt(rholay(i)*rholay(i-1))
         z = z + ((play(i-1)-play(i))/(glay(i)*drho))
         Rlay(i) = R0+z
      end do
      write(*,*) "Rlay calculated"
      do i=1,Nlayers
        zlay(i) = Rlay(i)-Rlay(Nlayers) 
        !print*,i,Rlay(i)/R0,zlay(i)
      enddo

*     --------------------------------------
*     ***  calculate Diffusion constant  ***
*     --------------------------------------
      write(*,*) "Calculating diffusion constant D ..."
      do i=Nlayers,1,-1
         Hplay = bk*Tlay(i)/(glay(i)*mulay(i)*amu)
         if (Tlay(i)>Teff) Hp=Hplay
         Difflay(i) = Hplay * ABS(Wmlay(i))
      enddo
      write(*,*) "diffusion constant D calculated"

      if (newMethod) then
*     --------------------
*     ***  new method  ***
*     --------------------
      !--- standard deviations of vertical velocities ---
      do i=1,Nlayers
        low = max(1,i-2)                 ! consider vz values on [i-2,i+3]
        up  = min(Nlayers,i+3)           ! -> sigma_vz between boxes i,i+1
        NNN = up-low+1
        meanv  = 0.d0
        meanv2 = 0.d0
        do j=low,up
          meanv  = meanv  + Wmlay(j) 
          meanv2 = meanv2 + (Wmlay(j))**2 
        enddo
        meanv  = meanv/NNN
        meanv2 = meanv2/NNN
        vz0(i) = SQRT(meanv2 - meanv**2) ! standard deviation
      enddo 
      !--- smooth these three times (box-car) ---
      do it=1,3
        do i=1,Nlayers
          low = max(1,i-1)
          up  = min(Nlayers,i+1)
          NNN = up-low+1
          meanv = 0.d0
          do j=low,up
            meanv = meanv + vz0(j)
          enddo
          vz(i) = meanv/NNN
          Hplay = bk*Tlay(i)/(glay(i)*mulay(i)*amu)
          Difflay(i) = Hplay * vz(i)
          if (it==3) print'(I4,0pF7.1,4(1pE11.3))',
     &       i,Tlay(i),play(i)/bar,Wmlay(i),vz(i),Difflay(i)
        enddo  
        vz0 = vz
      enddo  
      endif

*     -----------------------------
*     ***  add micro diffusion  ***     
*     -----------------------------
      !write(*,'(99(A10))') 'p[bar]','n[cm-3]','vth[km/s]',
     &!                     'l[cm]','Hp[cm]','Dmix','Dmicro'
      do i=Nlayers-1,1,-1
        Hplay  = bk*Tlay(i)/(glay(i)*mulay(i)*amu)
        if (Tlay(i)>Teff) Hp=Hplay
        ngas   = play(i)/bk/Tlay(i)
        lmean  = 1.d0/(2.1E-15*ngas)
        Kn     = lmean/Hp
        vth    = SQRT(8.d0*bk*Tlay(i)/(pi*mulay(i)*amu))
        Dmicro = 1.d0/3.d0*vth*lmean
        !write(*,'(99(1pE10.2))') play(i)/bar,ngas,vth/km,lmean,Hplay,
     >  !        Difflay(i),Dmicro
        Difflay(i) = Difflay(i) + Dmicro
      enddo  

      end
