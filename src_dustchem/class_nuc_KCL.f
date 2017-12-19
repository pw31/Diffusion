*********************************************************************
      SUBROUTINE CLASS_NUC_KCL(T,SS,nKCL,Jst,Nst)
*********************************************************************
*****                                                           *****
*****  berechnet die KCL-Keimbildungsrate                       *****
*****  mit Hilfe von klassischer Nukleationstheorie.            *****
*****                                                           *****
*****  INPUT:    T = Gastemperatur [K]                          *****
*****           SS = Uebersaettigungsverhaeltnis                *****
*****                supersaturation ratio                      *****
*****                                                           *****
*****  OUTPUT: Jst = Keimbildungsrate [cm^-3 s^-1]              *****
*****                nucleation rate                            *****
*****          Nst = Groesse des kritischen Clusters            *****
*****                critical cluster size                      *****
*****                                                           *****
*********************************************************************
      use NATURE,ONLY: amu,pi,bk
      implicit none
      real*8,intent(IN) :: T,SS,nKCL
      real*8,intent(OUT) :: Jst,Nst
      real*8 :: slog,thetun,thetaN,x0,x1,x2,x3,dgdn,fst
      real*8 :: zeldof,vth,beta,fNst
      logical :: IS_NAN

*     --------------------------------------------------------------
*     ***  Materialkonstanten von KCl(fest):                     ***
*     ***  sigma = Oberflaechenspannung in erg/cm^2              ***
*     ***        = 110 erg/cm^2 (Graham Lee)                     ***
*     ***  Nf = Clustergroesse, bei der sich Theta halbiert      ***
*     ***  f0 = 4*pi*a0^2 Monomeroberflaeche in cm^2             ***
*     ***  alfa = Sticking-coefficients                          ***
*     --------------------------------------------------------------
      real*8,parameter :: sigma=110.d0, Nf=0.d0
      real*8,parameter :: f0=7.59d-15, molg=74.5d0, alfa=1.d0

      slog = DLOG(SS)
      if (SS.le.1.d0) then
        Jst = 0.d+0
        Nst = 9.d+99
        return
      endif

*     ---------------------------------------------------------------
*     ***  Groesse des krit. Clusters nach dem Troepfchenmodel    ***
*     ***  Size of critical cluster according to droplet model    ***
*     ---------------------------------------------------------------
      thetun = f0*sigma/bk
      x0     = 2.d0*thetun/(3.d0*T*slog)
      x1     = Nf**(1.d0/3.d0)
      x2     = x1/x0
      x3     = 0.5d0*(1.d0+DSQRT(1.d0+2.d0*x2)) - x2
      Nst    = 1.d0 + (x0*x3)**3
      if (Nst.le.1.d0) Nst=1.000000001d0
*
*     ------------------------------------------------
*     ***  Teilchenhaeufigkeit des krit. Clusters  ***
*     ***  Number density of critical cluster      ***
*     ------------------------------------------------
      x0     = x0*slog
      x2     = (Nst-1.d0)**(1.d0/3.d0) + x1
      x3     = 1.d0/(Nst-1.d0)**(2.d0/3.d0)
      dgdn   = x0*x3* ( 1.d0/x2**2 + x1/x2**3 ) / 3.d0
      zeldof = DSQRT(dgdn/(2.d0*pi))
      thetaN = thetun/(1.d0+(Nf/(Nst-1.d0))**(1.d0/3.d0))
      x1     = (Nst-1.d0)*slog - (thetaN/T)
     &         *(Nst-1.d0)**(2.d0/3.d0)
      fst    = nKCL*DEXP(x1)

*     ----------------------------------
*     ***  Wachstumsgeschwindigkeit  ***
*     ***  Growth velocity           ***
*     ----------------------------------
      vth  = DSQRT(bk*T/(2.d0*pi*molg*amu))
      beta = vth*alfa*nKCL
      fNst = f0*Nst**(2.d0/3.d0)

*     --------------------------
*     ***  Keimbildungsrate  ***
*     ***  Nucleation rate   ***
*     --------------------------
      Jst = fst*beta*fNst*zeldof

#ifdef CHECK_NAN
      if (IS_NAN(Jst)) then
        print*,"*** NaN in class_nuc_KCl."
        print*,T,SS,slog
        print*,nKCl,vth
        print*,x0,x1,x2,x3
        print*,dgdn,zeldof,thetaN
        print*,f0
        print*,Nst,Jst
        stop
      endif   
#endif

      RETURN
      end
