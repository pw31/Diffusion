************************************************************************
      SUBROUTINE CREATE_STRUCTURE
************************************************************************
      use NATURE,ONLY: bar,bk,amu,km,pi,mel,grav,Mearth
      use PARAMETERS,ONLY: Tcrust,pmin,pmax,logg,Rplanet,Hp,
     >                     vzconst,pconst,beta
      use READMODEL,ONLY: Rlay,Tlay,play,rholay,glay,vconvlay,
     >                    zlay,mulay,Difflay,Nlayers
      use ELEMENTS,ONLY: NELEM,elnr,elcode,elnam,eps0,eps_solar,mass,muH
      use CHEMISTRY,ONLY: NELM,NMOLE,elnum,cmol,catm,el,charge,molmass
      use DUST_DATA,ONLY: NDUST
      use EXCHANGE,ONLY: nel,nat,nion,nmol,chi,inactive,N,O
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer :: i,it,iz
      real :: pwork,pp,nH,Tg,gg,Hplay,pconv,grad,ngas,lmean,Dmicro
      real :: sumn,sumnm,mu,Kn,vth,vz,TATMOS,Mpl,zz,dz
      real(kind=qp) :: eps(NELEM)
      integer,dimension(1000) :: flag_conv,Z
      logical :: conv
      character :: CR = CHAR(13)

      if (.not.allocated(nmol)) then
        allocate(nmol(NMOLE),chi(NDUST),inactive(NMOLE))
        inactive = .false.
      endif  

      print*
      print'(" create structure  Tcrust=",0pF6.1," K  psurf=",
     >       0pF7.3," bar")',Tcrust,pmax
      print*,"====================================================="

      !--- solve hydrostatic equilibrium ---
      eps_solar(N) = 80.0
      eps_solar(O) = 20.0
      eps0= eps_solar
      eps = eps_solar
      pp  = 2.0*pmax*bar
      mu  = 2.3*amu
      zz  = 0.d0
      gg  = 10.d0**logg
      Mpl = gg*Rplanet**2/grav
      muH = 28*amu
      !print'(A4,A12,A12,A10)',"iz","z[km]","p[bar]","T[K]"

      do iz=5000,1,-1
        Tg = TATMOS(pp)
        nH = mu*pp/(bk*Tg)/muH
        !--- iterate chemistry for given pressure ---
        do it=1,99
          call GGCHEM(nH,Tg,eps,.false.,0)
          sumn  = nel
          sumnm = nel*mel
          do i=1,NELEM
            if (nat(i)==0.Q0) cycle
            sumn  = sumn  + nat(i)
            sumnm = sumnm + nat(i)*mass(i)
          enddo  
          do i=1,NMOLE
            sumn  = sumn  + nmol(i)
            sumnm = sumnm + nmol(i)*molmass(i)
          enddo
          muH   = sumnm/nH
          mu    = sumnm/sumn
          pwork = sumn*bk*Tg
          nH    = nH*pp/pwork
          if (ABS(pp/pwork-1.d0)<1.E-8) exit
        enddo  
        write(*,'(".",$)') 
        !print'(I4,2(1pE12.3),0pF10.2)',iz,zz/km,pp/bar,Tg
        glay(iz) = gg
        zlay(iz) = zz
        Rlay(iz) = Rplanet+zz
        zlay(iz) = zz
        Tlay(iz) = Tg
        play(iz) = pp
        rholay(iz) = sumnm
        mulay(iz) = mu/amu
        Nlayers = 5001-iz
        if (pp<pmin*bar) exit      
        
        gg = grav*Mpl/Rlay(iz)**2
        Hp = bk*Tg/(gg*mu) 
        if (iz==1) then
          print'("  planet mass =",0pF10.2," Mearth")',Mpl/Mearth
          print'(" scale height =",0pF10.2," km")',Hp/km
        endif  
        dz = Hp/50.0
        zz = zz+dz
        pp = EXP(LOG(pp)-dz/Hp)
        Tg = TATMOS(pp)
        nH = mu*pwork/(bk*Tg)
      enddo
      print*
      print'(I4," layers created. pmin/bar =",1pE11.4," pmax/bar =",
     >       1pE11.4)',Nlayers,pp/bar,2*pmax
      print'(4x,"atmosphere extension dR/R =",1pE11.4)',zz/Rplanet

*     -------------------------
*     ***  shift to bottom  ***
*     -------------------------
      do i=1,Nlayers
        iz = 5000-Nlayers+i
        glay(i) = glay(iz)
        zlay(i) = zlay(iz)
        Rlay(i) = Rlay(iz)
        zlay(i) = zlay(iz)
        Tlay(i) = Tlay(iz)
        play(i) = play(iz)
        rholay(i) = rholay(iz)
        mulay(i) = mulay(iz)
        !print*,i,Rlay(i),zlay(i)
      enddo  

*     -----------------------------------------
*     ***  calculate Diffusion coefficient  ***
*     -----------------------------------------
      do i=1,Nlayers
        vz    = vzconst * MIN(1.d0,play(i)/pconst)**beta
        Hplay = bk*Tlay(i)/(glay(i)*mulay(i)*amu) 
        Difflay(i) = vz*Hplay
      enddo

*     -----------------------------
*     ***  add micro diffusion  ***     
*     -----------------------------
      !write(*,'(99(A10))') 'p[bar]','n[cm-3]','vth[km/s]',
     &!                     'l[cm]','Hp[cm]','Dmix','Dmicro'
      do i=Nlayers-1,1,-1
        Hplay  = bk*Tlay(i)/(glay(i)*mulay(i)*amu)
        ngas   = play(i)/bk/Tlay(i)
        lmean  = 1.d0/(2.1E-15*ngas)
        Kn     = lmean/Hp
        vth    = SQRT(8.d0*bk*Tlay(i)/(pi*mulay(i)*amu))
        Dmicro = 1.d0/3.d0*vth*lmean
        !write(*,'(99(1pE10.2))') play(i)/bar,ngas,vth/km,lmean,Hplay,
     >  !        Difflay(i),Dmicro
        Difflay(i) = Difflay(i) + Dmicro
      enddo  

    3 format(2X, I2, 2X, 1D10.3, 2x, f10.4, 3(2X, 1D10.3), 2x, f6.3,
     >       3(2X, 1PE10.3))
    4 format(3X, I2, 13(1X, 1D8.3), 1x, I2)
  100 format(f12.3, 1x, f12.3, 1x, f12.3)
  200 format(i5)
  210 format(8(i5, 1x))
  300 format(8(e15.8, 1x)) 
 1000 format(' Teff=',0pF8.3,' logg=',0pF5.2,' mixLengthPara=',0pF5.2)
 1100 format(i4,0pF11.3,1pE11.3,1x,0pF11.2,1x,2(1pe11.3,1x)) 
 1110 format(a4,5(a11,1x)) 
 1200 format(i4,"   Hp[km]=",F8.3,"   vconv[m/s]=",F8.3)
      end


***********************************************************************
      real*8 function TATMOS(pdyn)
***********************************************************************
***  parameteric T(p) structure according to Blecic+2017            ***
***  (http://adsabs.harvard.edu/abs/2017ApJ...848..127B),           ***
***  Appendix A, Parametrization Scheme II. In her paper,           ***
***  primary parameters are p0,p1,p2,p3,T3,a1,a2, whereas here      ***
***  primary parameters are p0,p1,p2,p3,T0,T1,T3.                   ***  
***********************************************************************
      use NATURE,ONLY: bar
      use PARAMETERS,ONLY: pmin,pmax,Tcrust  ! [bar],[bar],[K]
      implicit none
      real*8,intent(IN) :: pdyn              ! p[dyn/cm2]
      real*8 :: p,dp,Tsum
      real*8,save :: T0,T1,T2,T3,p0,p1,p2,p3,a1,a2 
      logical,save :: firstCall=.true.
      integer :: j,nsmooth

      if (firstCall==.true.) then
        p3 = pmax      ! surface 
        p2 = 0.0       ! tropopause - computed below
        p1 = p3*1.E-5  ! stratopause
        p0 = pmin      ! exobase

        T3 = Tcrust    ! surface
        T2 = 0.7*T3    ! tropopause
        T1 = 2.0*T3    ! stratopause      
        T0 = T1-0.1    ! exobase

        a1 = LOG(p1/p0)/SQRT(T1-T0)
        a2 = LOG(p3/p1)/(SQRT(T1-T2)+SQRT(T3-T2))
        p2 = p1*EXP(a2*SQRT(T1-T2))
        firstCall = .false.
      endif  

      dp = LOG(p3/p0)/1000.0
      Tsum = 0.0
      nsmooth = 50
      do j=-nsmooth,nsmooth
        p = EXP(LOG(pdyn/bar)+j*dp)
        if (p<p1) then
          Tsum = Tsum + T0 + (LOG(p/p0)/a1)**2
        else if (p<p3) then
          Tsum = Tsum + T2 + (LOG(p/p2)/a2)**2
        else
          Tsum = Tsum + T3
        endif  
      enddo

      TATMOS = Tsum/(2.0*nsmooth+1.0)
      end
