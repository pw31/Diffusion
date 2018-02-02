*********************************************************************
      SUBROUTINE CLASS_NUC_SIO(T,SS,nSiO,Jst,Nst)
*********************************************************************
*****                                                           *****
*****  berechnet die SiO-Keimbildungsrate                       *****
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
      implicit none
      real*8,intent(in) :: T,SS,nSiO
      real*8,intent(out) :: Jst,Nst
      real*8 :: slog,slog2,T3

*     -----------------------------------------------
*     ***  for more details see class_nuc_TiO2.f  ***
*     -----------------------------------------------
      !print*,"T,SS,nSiO=",T,SS,nSiO
      slog = DLOG(SS)
      if (SS.le.1.d0) then
        Jst = 0.d+0
        Nst = 9.d+99
        return
      end if  

*     -----------------------------------
*     ***  Keimbildungsrate           ***
*     ***  Nucleation rate            ***
*     ***  Eq 16 in Gail et al. 2013  ***
*     -----------------------------------
      Nst   = 1.d0                   ! not calculated here
      T3    = T*T*T
      slog2 = slog*slog  
      Jst   = nSiO*nSiO * DEXP(1.33d0 - 4.4d12/(T3*slog2))

      end
