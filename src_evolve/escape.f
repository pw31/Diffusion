************************************************************************
      subroutine ESCAPE
************************************************************************
      use PARAMETERS,ONLY: implicit,Rplanet,logg,verbose
      use GRID,ONLY: zz,Npoints
      use STRUCT,ONLY: Temp,nHtot,nHeps
      use ELEMENTS,ONLY: NELEM,elnam,eps0,mass
      use CHEMISTRY,ONLY: NMOLE,molmass,elnum,cmol,el,elnum,NELM,
     >                    m_kind,m_anz
      use EXCHANGE, ONLY: nat,nmol,nel
      use NATURE,ONLY: bk,pi
      use JEANS_ESCAPE, ONLY: jpern

      implicit none

      real*8 :: nH,Tg,ztop,vesc,gravity,vth,JJJ
      integer :: i,j,e
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: eps(NELEM),flux(NELEM)
      real(kind=qp) :: vmol_top(NMOLE),Jeans_mol(NMOLE),Jeans_at(NELEM)
&,vat_top(NELEM)

      !--------------------------------------------------------------------------
      ! ***  call GGchem on top of atmosphere to get molecular particle dens. ***
      !--------------------------------------------------------------------------
      nH  = nHtot(Npoints)               ! density at top of atmosphere
      Tg  = Temp(Npoints)                ! temperature at top of atmosphere
      eps = eps0                         ! set element abundances to default
      do i=1,NELEM                    
        eps(i) = nHeps(i,Npoints)/nH     ! use element abundance at top
      enddo         
      call GGCHEM(nH,Tg,eps,.false.,0) 
      
      !--------------------------------------------------------------------------
      ! ***  Escape velocity of the planet at top of atmosphere. ***
      !--------------------------------------------------------------------------

      ztop = zz(Npoints) !Height of top of atmosphere in meters [cm]
      vesc = Sqrt((2.0*10**logg*Rplanet**2)/(Rplanet+ztop)) !Escape vel of planet

      !--------------------------------------------------------------------------
      ! ***  Jeans Escape ***
      !--------------------------------------------------------------------------

      flux=0.0
      JJJ = 0.0

      do i=1,NMOLE !Over all molecules
        vth = Sqrt((2.0*bk*Tg)/molmass(i)) !molecule vel at top of atmos
        JJJ = ((nmol(i)*vth)/(2.0*Sqrt(pi)))*((vesc/vth)**2+1.0)*
&Exp(-(vesc/vth)**2)
        do j=1,m_kind(0,i)
          e = m_kind(j,i) 
          el = elnum(e)
          flux(el) = flux(el) + m_anz(j,i)*JJJ 
          print*,cmol(i),el,m_anz(j,i),JJJ
      	enddo
      enddo

      do el=1,NELEM !Over all atoms
      	vth  = Sqrt((2.0*bk*Tg)/mass(i)) !atom vel at top of atmos
	JJJ = ((nat(i)*vth)/(2.0*Sqrt(pi)))*((vesc/vth)**2+1.0)*
&Exp(-(vesc/vth)**2)
        do j=1,m_kind(0,i)
          e = m_kind(j,i) 
          flux(el) = flux(el) + JJJ 
          print*,cmol(i),el
      	enddo
      enddo

      do el=1,NELEM
        jpern(el) = flux(el)/(nH*eps(el))
      enddo

      print*,"jpern",jpern
      stop

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

