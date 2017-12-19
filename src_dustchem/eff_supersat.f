*********************************************************************
      SUBROUTINE EFF_SUPERSAT(nH,T,eps,bmix,Sat,effSat)
*********************************************************************
*****                                                           *****
*****  EINGABE:  nH = total hydrogen density [cm^-3]            *****
*****             T = Gastemperatur [K]                         *****
*****       bmix(:) = Vs/Vtot                                   *****
*****                                                           *****
*****  AUSGABE: effSat = effektive Uebersaettigungsverhaelt.    *****
*****                                                           *****
*********************************************************************
      use NATURE,ONLY: pi,bk
      use ELEMENTS,ONLY:  NELEM,qp
      use CHEMISTRY,ONLY: NMOLE
      use SPECIES,ONLY: NSPECIES,spnam,spnr,spmass,keysp
      use REACTIONS,ONLY: NREAC,neduct,nprod,reac_sp,reac_nu,stick
      use DUST_DATA,ONLY: NDUST,dust_vol,dust_nam
      use EXCHANGE,ONLY: nel,nat,nion,nmol,chi
      implicit none      
      real*8,intent(IN) :: nH,T,bmix(NDUST)
      real(kind=qp),intent(IN):: eps(NDUST)
      real(kind=qp),intent(OUT):: Sat(NDUST)
      real*8,intent(OUT):: effSat(NDUST)
      integer :: r,i,j,sp,dustnr
      real(kind=qp) :: sum1(NDUST),sum2(NDUST)
      real(kind=qp) :: nsp(NSPECIES),vrel(NSPECIES),Sr,sum,delV
      real(kind=qp) :: msum,rsum,eq_fak,rate,rat1,stoi,mr,dustnu

*     --------------------------------------
*     ***  calculate particle densities  ***
*     --------------------------------------
      call GGCHEM(nH,T,eps,.false.,0)
      do i=1,NSPECIES
        j = spnr(i)
        if (j.gt.1000) then
          nsp(i) = nat(j-1000)
        else
          nsp(i) = nmol(j)
        endif
        vrel(i) = SQRT(bk*T/(2.d0*pi*spmass(i)))
      enddo
*
*     ----------------------------------------------
*     ***  calculate pure supersaturation ratio  ***
*     ----------------------------------------------
      call SUPERSAT(T,nat,nmol,Sat)

*     --------------------------------
*     ***  growth and evaporation  ***
*     --------------------------------
      sum1 = 0.d0
      sum2 = 0.d0
      do r=1,NREAC
        dustnr = reac_sp(r,1,2)
        dustnu = reac_nu(r,1,2)       
        rate   = 9.d+99
        eq_fak = 9.d+99
        msum   = 0.d0
        rsum   = 0.d0
        do j=1,neduct(r)
          sp = reac_sp(r,j,1)                       ! index for growth species
          if (keysp(sp)) then
	    stoi = 1.d0/DBLE(reac_nu(r,j,1))        ! ratio of stoichiom. coeffs
	    rat1 = stick(r)*nsp(sp)*vrel(sp)*stoi   ! surface reac/cm^2/s
            rsum = rsum + 1.Q0/rat1**2
            msum = msum + (stoi*dustnu)/rat1**2
          endif
        enddo
        mr = msum/rsum
        Sr = Sat(dustnr)**mr                        ! reaction supersat.ratio
        eq_fak = bmix(dustnr)/Sr                    ! factor for evap-rate
        sum1(dustnr) = sum1(dustnr) + dustnu * rate * eq_fak
        sum2(dustnr) = sum2(dustnr) + dustnu * rate
      enddo
*
*     ----------------------------------------------------
*     ***  calculate effective supersaturation ratios  ***
*     ----------------------------------------------------
      do i=1,NDUST
        effSat(i) = sum2(i)/sum1(i)
c       if (sum2(i)-sum1(i).ne.0.d0) then 
c         effSat(i) = 1.d0/(1.d0-sum1(i)/sum2(i))
c         effSat(i) = sum2(i)/(sum2(i)-sum1(i))
c       else
c         effSat(i) = Sat(i)
c       endif  
        !print*,dust_nam(i),REAL(Sat(i)),effSat(i)
      enddo

      RETURN 
      END 
