************************************************************************
      subroutine ESCAPE(verbose)
************************************************************************
      use PARAMETERS,ONLY: implicit,Rplanet,logg
      use GRID,ONLY: zz,Npoints
      use STRUCT,ONLY: Temp,nHtot,nHeps
      use ELEMENTS,ONLY: NELEM,elnam,eps0,mass
      use CHEMISTRY,ONLY: NMOLE,molmass,elnum,cmol,elnum,NELM,
     >                    m_kind,m_anz,iel=>el
      use EXCHANGE, ONLY: nat,nmol,nel
      use NATURE,ONLY: bk,pi,grav,km
      use JEANS_ESCAPE, ONLY: jpern
      implicit none
      integer,intent(in) :: verbose
      real*8 :: nH,Tg,Mpl,ztop,vesc,gravity,vth,JJJ
      integer :: i,j,e,el
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
      
      !-----------------------------------------------
      ! ***  Escape velocity at top of atmosphere  ***
      !-----------------------------------------------
      ztop = zz(Npoints)                   ! height of top of atmosphere [cm]
      Mpl  = 10.0**logg*Rplanet**2/grav    ! mass of planet [g]
      vesc = Sqrt(2.0*grav*Mpl/(Rplanet+ztop)) 
      if (verbose>0) print*,"escape velocity[km/s] =",vesc/km

      !-----------------------
      ! ***  Jeans Escape  ***
      !-----------------------
      flux(:) = 0.0                        ! sum_k j_k*stoich_k,el
      do i=1,NMOLE                         ! loop over molecules
        vth = SQRT((2.0*bk*Tg)/molmass(i)) ! thermal velocity
        JJJ = nmol(i)*vth/(2.0*Sqrt(pi))*((vesc/vth)**2+1.0)
     >        * EXP(-(vesc/vth)**2)        ! Feng+2015, Eq(1)
        if (verbose>0) print'(A12,3(1pE11.3))',
     >                 cmol(i),nmol(i),vth/km,JJJ/nmol(i)
        do j=1,m_kind(0,i)                 ! loop of elements in molecule
          e = m_kind(j,i)                  ! index in ggchem-list
          if (e==iel) cycle
          el = elnum(e)                    ! index of element
          flux(el) = flux(el) + m_anz(j,i)*JJJ 
          !print*,cmol(i),elnam(el),m_anz(j,i),JJJ
      	enddo
      enddo
      do e=1,NELM                          ! loop over atoms
        if (e==iel) cycle                  ! scip if electrons
        el = elnum(e)
      	vth = Sqrt((2.0*bk*Tg)/mass(el))   ! thermal velocity
	JJJ = nat(el)*vth/(2.0*SQRT(pi))*((vesc/vth)**2+1.0)
     >        * EXP(-(vesc/vth)**2)
        flux(el) = flux(el) + JJJ 
        if (verbose>0) print'(A2,3(1pE11.3))',
     >                 elnam(el),nat(el),vth/km,JJJ/nat(el)
      enddo

      !-------------------------------------------------------------
      ! ***  compute escaping flux per element particle density  ***
      !-------------------------------------------------------------
      do e=1,NELM
        if (e==iel) cycle                  ! scip if electrons
        el = elnum(e) 
        jpern(el) = flux(el)/(nH*eps(el))
        print*,elnam(el)//" jpern=",jpern(el)
      enddo
      if (verbose>0) stop

      !--------------------------------------------------------------------------
      ! ***  Jeans Escape plot ***
      !--------------------------------------------------------------------------
      !Writing out Jeans escape per molecule for plotting
      !open(unit=75,file='jeans.dat',status='replace')
      !do i=1,NMOLE
      !	write(75,*) (trim(cmol(i))),Jeans_mol(i),Jeans_at(i)
      !enddo
      !close(75)

      end subroutine ESCAPE

