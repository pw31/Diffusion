************************************************************************
      subroutine ESCAPE(deltat,reduced,verbose)
************************************************************************
      use PARAMETERS,ONLY: implicit,Rplanet,logg
      use GRID,ONLY: zz,Npoints,zweight
      use STRUCT,ONLY: Temp,nHtot,nHeps
      use ELEMENTS,ONLY: NELEM,elnam,eps0,mass
      use CHEMISTRY,ONLY: NMOLE,molmass,elnum,cmol,elnum,NELM,
     >                    m_kind,m_anz,iel=>el
      use EXCHANGE, ONLY: nat,nmol,nel,H
      use NATURE,ONLY: bk,pi,grav,km,yr
      use JEANS_ESCAPE,ONLY: Ttop,Hfrac,jpern
      implicit none
      real*8,intent(inout) :: deltat
      logical,intent(out) :: reduced
      integer,intent(in) :: verbose
      real*8 :: nH,Tg,Mpl,ztop,vesc,gravity,vth,JJJ,Natmos(NELEM)
      real*8 :: tau,Hsum
      integer :: i,j,e,el,STINDEX,iH2,iH2O
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: eps(NELEM),flux(NELEM)

      !---------------------------------------------------------------------
      ! ***  call GGchem on top of atmosphere to get particle densities  ***
      !---------------------------------------------------------------------
      nH  = nHtot(Npoints)                 ! density at top of atmosphere
      Tg  = Temp(Npoints)                  ! temperature at top of atmosphere
      eps = eps0                           ! set element abundances to default
      do i=1,NELEM                    
        eps(i) = nHeps(i,Npoints)/nH       ! use element abundance at top
      enddo         
      call GGCHEM(nH,Tg,eps,.false.,0) 
      iH2  = STINDEX(cmol,NMOLE,'H2')
      iH2O = STINDEX(cmol,NMOLE,'H2O')
      Hfrac(1) = 1.0*nat(H)
      Hfrac(2) = 2.0*nmol(iH2)
      Hfrac(3) = nH*eps(H)-Hfrac(1)-Hfrac(2) ! all other H-species escape like H2O
      Hsum  = nH*eps(H)
      Hfrac = Hfrac/Hsum
      if (verbose>0) print'(" Hfrac=",3(1pE12.4))',Hfrac
      
      !-----------------------------------------------
      ! ***  Escape velocity at top of atmosphere  ***
      !-----------------------------------------------
      Tg   = Ttop                          ! artificially high T 
      ztop = zz(Npoints)                   ! height of top of atmosphere [cm]
      Mpl  = 10.0**logg*Rplanet**2/grav    ! mass of planet [g]
      vesc = Sqrt(2.0*grav*Mpl/(Rplanet+ztop)) 
      if (verbose>1) print*,"escape velocity[km/s] =",vesc/km

      !-----------------------
      ! ***  Jeans Escape  ***
      !-----------------------
      flux(:) = 0.0                        ! sum_k j_k*stoich_k,el
      do e=1,NELM                          ! loop over atoms
        if (e==iel) cycle                  ! skip if electrons
        el = elnum(e)                      ! index of element
      	vth = Sqrt((2.0*bk*Tg)/mass(el))   ! thermal velocity
	JJJ = nat(el)*vth/(2.0*SQRT(pi))*((vesc/vth)**2+1.0)
     >        * EXP(-(vesc/vth)**2)        ! Feng+2015, Eq(1)  [1/cm2/s]
        flux(el) = flux(el) + JJJ 
        if (el==H) jpern(NELEM+1)=JJJ/nat(el)
        if (verbose>1) print'(A2,4(1pE11.3))',elnam(el),
     >                 nat(el),vth/km,JJJ/nat(el),JJJ/(nH*eps(el))
      enddo
      do i=1,NMOLE                         ! loop over molecules
        vth = SQRT((2.0*bk*Tg)/molmass(i)) ! thermal velocity
        JJJ = nmol(i)*vth/(2.0*Sqrt(pi))*((vesc/vth)**2+1.0)
     >        * EXP(-(vesc/vth)**2)        ! Feng+2015, Eq(1)  [1/cm2/s]
        if (i==iH2)  jpern(NELEM+2)=JJJ/nmol(i)
        if (i==iH2O) jpern(NELEM+3)=JJJ/nmol(i)
        do j=1,m_kind(0,i)                 ! loop over elements in molecule
          e = m_kind(j,i)                  ! index in ggchem-list
          if (e==iel) cycle                ! skip if electrons
          el = elnum(e)                    ! index of element
          flux(el) = flux(el) + m_anz(j,i)*JJJ 
          if (verbose>1) print'(A12,4(1pE11.3))',cmol(i),
     >       nmol(i),vth/km,JJJ/nmol(i),m_anz(j,i)*JJJ/(nH*eps(el))
      	enddo
      enddo

      !-------------------------------------------------------------
      ! ***  compute escaping flux per element particle density  ***
      !-------------------------------------------------------------
      if (verbose>0) print'(" H jpern contributions:",3(1pE12.3))',
     >       jpern(NELEM+1)*nat(H)/(nH*eps(H)),
     >       jpern(NELEM+2)*nmol(iH2)*2/(nH*eps(H)),
     >       jpern(NELEM+3)*nmol(iH2O)*2/(nH*eps(H))
      tau = 9.E+99
      do e=1,NELM
        if (e==iel) cycle                  
        el = elnum(e) 
        jpern(el) = flux(el)/(nH*eps(el))                   ! [cm/s]
        Natmos(el) = 0.d0
        do i=0,Npoints
          Natmos(el) = Natmos(el) + nHeps(el,i)*zweight(i)  ! [1/cm2]
        enddo
        tau = MIN(tau,Natmos(el)/flux(el))                  ! [s]
        if (verbose>1) print'(A3,"  jpern=",1pE10.3,"  tau[yrs]=",
     >       1pE10.3)',elnam(el),jpern(el),Natmos(el)/flux(el)/yr
      enddo

      !---------------------------
      ! ***  timestep control  ***
      !---------------------------
      reduced = .false.
      !if (deltat>2.0*tau) then
      !  deltat = 2.0*tau
      !  reduced = .true.
      !  print*,"*** ESCAPE: timestep too large"
      !endif  
      if (verbose>0) then
        print'("ESCAPE: dt,tau[yrs]=",2(1pE11.3))',deltat/yr,tau/yr
        print'("expected decrease of NHatmos:  dt=",1pE10.3,"  NH=",
     >       1pE18.11," ->",1pE18.11)',deltat,
     >       Natmos(H),Natmos(H)-flux(H)*deltat
      endif  

      end subroutine ESCAPE

