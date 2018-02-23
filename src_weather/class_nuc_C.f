*********************************************************************
      SUBROUTINE CLASS_NUC_C(T,SS,nC1,nC2,nC2H,nC2H2,nC3,Jst,Nst)
*********************************************************************
*****                                                           *****
*****  berechnet die Kohlenstoff-Keimbildungsrate               *****
*****  mit Hilfe von klassischer Nukleationstheorie.            *****
*****                                                           *****
*****  EINGABE:   T = Gastemperatur [K]                         *****
*****            SS = Uebersaettigungsverhaeltnis               *****
*****           nC1 = Teilchendichte des Monomers C1 [cm^-3]    *****
*****           nC2 = Teilchendichte des Dimers   C2 [cm^-3]    *****
*****           nC3 = Teilchendichte des Trimers  C3 [cm^-3]    *****
*****                                                           *****
*****  AUSGABE: Jst = Keimbildungsrate [cm^-3 s^-1]             *****
*****           Nst = Groesse des kritischen Clusters           *****
*****                                                           *****
*********************************************************************
       use NATURE,ONLY: pi,bk,amu
       implicit none
       real*8,intent(in) :: T,SS,nC1,nC2,nC2H,nC2H2,nC3
       real*8,intent(out) :: Jst,Nst
       real*8,parameter :: mH=1.008, mC=12.01
       real*8,parameter :: mC1=mC, mC2=2*mC, mC3=3*mC
       real*8,parameter :: mC2H=2*mC+mH, mC2H2=2*mC+2*mH
       real*8 :: Tvar,slog,thetun,thetaN,x0,x1,x2,x3,dgdn,fst
       real*8 :: zeldof,vth,beta,Ast

       !--------------------------------------------------------------
       !***  Materialkonstanten von Graphit:                       ***
       !***  sigma = Oberflaechenspannung in [erg/cm2]             ***
       !***  Nl = Clustergroesse, bei der sich Theta halbiert      ***
       !***  f0 = 4*pi*a0^2 Monomeroberflaeche in [cm2]            ***
       !***  alf1,alf2,alf3 = Sticking-coefficients fuer C1,C2,C3  ***
       !***    ((((( nach GAIL:   sigma=1400   Nl=5 )))))          ***
       !-------------------------------------------------------------- 
       real*8,parameter :: sigma=1400.d0, Nl=5.d0, f0=2.07e-15
       real*8,parameter :: alf1=0.37d0, alf2=0.34d0, alf3=0.08d0

       slog = DLOG(SS)
       if (SS.le.1.d0) then
         Jst = 0.d0
         Nst = 9.d+99
         return
       endif

       !--------------------------------------------------------------
       !***  Groesse des krit. Clusters nach dem Troepfchenmodel   ***
       !--------------------------------------------------------------
       thetun = f0*sigma/bk
       x0     = 2.d0*thetun/(3.d0*Tvar*slog)
       x1     = Nl**(1.d0/3.d0)
       x2     = x1/x0
       x3     = 0.5d0*(1.d0+SQRT(1.d0+2.d0*x2)) - x2
       Nst    = 1.d0 + (x0*x3)**3 
       if (Nst.le.1.d0) Nst=1.000000001d0

       !------------------------------------------------
       !***  Teilchenhaeufigkeit des krit. Clusters  ***
       !------------------------------------------------
       x0     = x0*slog
       x2     = (Nst-1.d0)**(1.d0/3.d0) + x1
       x3     = 1.d0/(Nst-1.d0)**(2.d0/3.d0)
       dgdn   = x0*x3* ( 1.d0/x2**2 + x1/x2**3 ) / 3.d0
       zeldof = DSQRT(dgdn/(2.d0*pi))
       thetaN = thetun/(1.d0+(Nl/(Nst-1.d0))**(1.d0/3.d0))
       x1     = (Nst-1.d0)*slog - (thetaN/Tvar)*(Nst-1.d0)**(2.d0/3.d0)
       if (x1.lt.-300.0) then
         Jst = 0.d0
         Nst = 0.d0
         return
       endif
       fst = nc1*DEXP(x1)

       !----------------------------------
       !***  Wachstumsgeschwindigkeit  ***
       !----------------------------------
       vth  = DSQRT(bk*T/(2.d0*pi*amu))
       beta = vth*(alf1*nC1  /DSQRT(mC1)            
     &       +2.d0*alf2*nC2  /DSQRT(mC2)            
     &       +2.d0*alf2*nC2H /DSQRT(mC2H)           
     &       +2.d0*alf2*nC2H2/DSQRT(mC2H2)         
     &       +3.d0*alf3*nC3  /DSQRT(mC3) )

       !--------------------------
       !***  Keimbildungsrate  ***
       !--------------------------
       Ast = f0*Nst**(2.d0/3.d0)
       Jst = fst*beta*Ast*zeldof

       end
