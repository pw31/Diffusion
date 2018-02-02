*********************************************************************
      SUBROUTINE CLASS_NUC_TIO2(T,SS,nTiO2,Jst,Nst)
*********************************************************************
*****                                                           *****
*****  berechnet die TiO2-Keimbildungsrate                      *****
*****  mit Hilfe von klassischer Nukleationstheorie.            *****
*****                                                           *****
*****  INPUT:       T = Gastemperatur [K]                       *****
*****              SS = Uebersaettigungsverhaeltnis             *****
*****                   supersaturation ratio                   *****
*****                                                           *****
*****  OUTPUT: Jstern = Keimbildungsrate [cm^-3 s^-1]           *****
*****                   nucleation rate                         *****
*****          Nstern = Groesse des kritischen Clusters         *****
*****                   critical cluster size                   *****
*****                                                           *****
*********************************************************************
      use NATURE,ONLY: pi,bk,amu
      implicit none
      real*8,intent(in)  :: T,SS,nTiO2
      real*8,intent(out) :: Jst,Nst
      real*8 :: slog,thetun,thetaN,x0,x1,x2,x3,dgdn,fst
      real*8 :: zeldof,vth,beta,fNst
*
*     --------------------------------------------------------------
*     ***  Materialkonstanten von TiO2(fest):                    ***
*     ***  sigma = Oberflaechenspannung in erg/cm^2              ***
*     ***        = 620   erg/cm^2 Jeong                          ***
*     ***        = 480.6 erg/cm^2 Graham Lee                     ***
*     ***  Nf = Clustergroesse, bei der sich Theta halbiert      ***
*     ***  f0 = 4*pi*a0^2 Monomeroberflaeche in cm^2             ***
*     ***  alfa = Sticking-coefficients                          ***
*     --------------------------------------------------------------
      real*8,parameter :: sigma=480.d0, Nf=0.d0
      real*8,parameter :: f0=4.808362d-15, molg=79.898d0, alfa=1.d0

*     ---------------------------------------------------
*     ***  Dampfdruck von TiO2 (Gas) ueber TiO2(fest) ***
*     ***  nach JANAF-Tafeln, elektr. Version 1985    ***
*     ***  Fit zwischen 500K und 2500K  (Woitke)      ***
*     ---------------------------------------------------
      !write(*,*) 'nTiO2 = ',nTiO2,' Psat = ',psat,' SS = ',SS
      slog = DLOG(SS)
      if (SS.le.1.d0) then
        Jst = 0.d+0
        Nst = 9.d+99
        return
      end if  

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
      if (Nst.le.1.d0) Nst=1.000000001D0

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
      x1     = (Nst-1.d0)*slog - (thetaN/T)*(Nst-1.d0)**(2.d0/3.d0)
      fst    = nTiO2*DEXP(x1)

*     ----------------------------------
*     ***  Wachstumsgeschwindigkeit  ***
*     ***  Growth velocity           ***
*     ----------------------------------
      vth  = DSQRT(bk*T/(2.d0*pi*molg*amu))
      beta = vth*alfa*nTiO2
      fNst = f0*Nst**(2.d0/3.d0)

*     --------------------------
*     ***  Keimbildungsrate  ***
*     ***  Nucleation rate   ***
*     --------------------------
      Jst = fst*beta*fNst*zeldof

 500  continue
      RETURN
      end
