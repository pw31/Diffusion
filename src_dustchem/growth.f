*********************************************************************
      SUBROUTINE GROWTH(T,Sat,bmix,nat,nmol,chinet,chi,cr)
*********************************************************************
      use NATURE,ONLY: bk,pi,km
      use ELEMENTS,ONLY:  NELEM,qp
      use CHEMISTRY,ONLY: NMOLE
      use SPECIES,ONLY: NSPECIES,spnam,spnr,spmass,keysp
      use REACTIONS,ONLY: NREAC,neduct,nprod,reac_sp,reac_nu,stick
      use DUST_DATA,ONLY: NDUST,dust_vol,dust_nam
      use EXCHANGE,ONLY: ipoint
      implicit none
      real*8,intent(in) :: T
      real(kind=qp),intent(in) :: Sat(NDUST)
      real(kind=qp),intent(in) :: nat(NELEM),nmol(NMOLE)
      real*8,intent(in)  :: bmix(NDUST)
      real*8,intent(out) :: chinet,chi(NDUST),cr(NREAC)
      real*8 :: nsp(NSPECIES),vrel(NSPECIES)
      real*8 :: dustnu,rat1,rate,eq_fak,msum,rsum,stoi,mr,Sr,dchi
      integer :: sp,j,r,dustnr,key
      logical :: out=.false.
      real*8,save :: Afak
      logical,save :: firstCall=.true.
!$omp threadprivate(firstCall,Afak)

      if (firstCall) then
        Afak = (36.d0*pi)**(1.d0/3.d0)
        firstCall = .false.
      endif

      !--------------------------------------------------------------
      !***  fetch particle densities, compute thermal velocities  ***
      !--------------------------------------------------------------
      do sp=1,NSPECIES
        j = spnr(sp)
        if (j.gt.1000) then
          nsp(sp) = nat(j-1000)
        else
          nsp(sp) = nmol(j)
        endif
        vrel(sp) = DSQRT(bk*T/(2.d0*pi*spmass(sp)))
      enddo

      chinet = 0.d0
      chi(:) = 0.d0
      do r=1,NREAC
        dustnr = reac_sp(r,1,2)                     ! first product is dust species
        dustnu = reac_nu(r,1,2)                     ! no of solid units produced
        key    = 0
        rate   = 9.d+99
        msum   = 0.d0
        rsum   = 0.d0
        !--------------------------------------
        ! ***  identify key species and mr  ***
        !--------------------------------------
        do j=1,neduct(r)               
          sp = reac_sp(r,j,1)                       ! index of growth species
          if (keysp(sp)) then
	    stoi = 1.d0/DBLE(reac_nu(r,j,1))        ! ratio of stoichiom. coeffs
	    rat1 = stick(r)*nsp(sp)*vrel(sp)*stoi   ! surface reac/cm^2/s
	    if (rat1.lt.rate) then
	      rate = rat1                           ! minimum is key 
              key  = sp
            endif  
            rsum = rsum + 1.d0/rat1**2
            msum = msum + (stoi*reac_nu(r,1,2))/rat1**2
          endif
        enddo
        mr = msum/rsum
        !---------------------------------------------------------
        ! ***  reaction upersaturation ratio and book-keeping  ***
        !---------------------------------------------------------
        Sr = Sat(dustnr)**mr                        ! reaction supersat.ratio
        eq_fak = 1.d0 - bmix(dustnr)/Sr             ! equilibrium factor 
        cr(r)  = Afak * dustnu * rate * eq_fak      ! net rate [cm-2 s-1]
        dchi   = cr(r) * dust_vol(dustnr)           ! volume increase [cm/s]
        chinet = chinet + dchi
        chi(dustnr) = chi(dustnr) + dchi

        if (out) then
!$omp critical(out90)
          if (neduct(r).eq.1) then
            write(90,2011) ipoint,r,stick(r),
     &      (reac_nu(r,j,1),spnam(reac_sp(r,j,1)),j=1,neduct(r)),
     &      reac_nu(r,1,2), dust_nam(dustnr),
     &      (reac_nu(r,j,2),spnam(reac_sp(r,j,2)),j=2,nprod(r))
          else if (neduct(r).eq.2) then
            write(90,2021) ipoint,r,stick(r),
     &      (reac_nu(r,j,1),spnam(reac_sp(r,j,1)),j=1,neduct(r)),
     &      reac_nu(r,1,2), dust_nam(dustnr),
     &      (reac_nu(r,j,2),spnam(reac_sp(r,j,2)),j=2,nprod(r))
          else if (neduct(r).eq.3) then
            write(90,2031) ipoint,r,stick(r),
     &      (reac_nu(r,j,1),spnam(reac_sp(r,j,1)),j=1,neduct(r)),
     &      reac_nu(r,1,2), dust_nam(dustnr),
     &      (reac_nu(r,j,2),spnam(reac_sp(r,j,2)),j=2,nprod(r))
          endif  
          write(90,*) " key species = "//trim(spnam(key))
          write(90,'("  n=",1pE8.2,"  vth/km=",0pF5.3,"  rate=",1pE8.2,
     &           "  mr=",0pF4.2,"  Sr=",1pE8.2,"  bmix=",1pE8.2,
     &           "  eq_fak=",1pE9.2)')
     &         nsp(key),vrel(key)/km,rate,mr,Sr,bmix(dustnr),eq_fak
          write(90,*)
!$omp end critical(out90)
        endif  
      enddo

 2011 format(I3,I3,0pF5.2,1x,1(I2,1x,a8),22x,'->',I2,1x,a10,9(I2,1x,a8))
 2021 format(I3,I3,0pF5.2,1x,2(I2,1x,a8),11x,'->',I2,1x,a10,9(I2,1x,a8))
 2031 format(I3,I3,0pF5.2,1x,3(I2,1x,a8)    ,'->',I2,1x,a10,9(I2,1x,a8))
      end 
