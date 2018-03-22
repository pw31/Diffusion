************************************************************************
      subroutine ESCAPE
************************************************************************
      use PARAMETERS,ONLY: implicit,Rplanet,logg,verbose
      use GRID,ONLY: zz,Npoints
      use STRUCT,ONLY: Temp,nHtot,nHeps
      use ELEMENTS,ONLY: NELEM,elnam,eps0
      use CHEMISTRY,ONLY: NMOLE,molmass,elnum,cmol
      use EXCHANGE, ONLY: nat,nmol,nel
      use NATURE,ONLY: bk,pi
      use JEANS_ESCAPE, ONLY: NMOLE,nel_top,nat_top,nmol_top

      implicit none

      real*8 :: nH,Tg,ztop,vesc,gravity
      integer :: i,j
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: eps(NELEM)
      real(kind=qp) :: vmol_top(NMOLE),Jeans(NMOLE)

      call READ_PARAMETER !Reading in mass of particles 

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
      nmol_top = nmol                    ! store results in module JEANS_ESCAPE
      nat_top  = nat
      nel_top  = nel

      
      !--------------------------------------------------------------------------
      ! ***  Escape velocity of the planet at top of atmosphere. ***
      !--------------------------------------------------------------------------

      !call READ_PARAMETER !Using values from input file
      ztop = zz(Npoints) !Height of top of atmosphere in meters [cm]
      vesc = Sqrt((2.0*10**logg*(Rplanet*10**5.0)**2.0)
&/(Rplanet*10**5.0+ztop)) !Escape vel of planet

      !--------------------------------------------------------------------------
      ! ***  Particle velocity at top of atmosphere. ***
      !--------------------------------------------------------------------------

      do i=1,NMOLE !Over all molecules
	vmol_top(i) = Sqrt((2.0*bk*Tg)/molmass(i)) !Particle vel at top of atmos
      enddo

      !--------------------------------------------------------------------------
      ! ***  Jeans Escape ***
      !--------------------------------------------------------------------------

      do i=1,NMOLE !Sum over all molecules
	Jeans(i) = ((nmol(i)*vmol_top(i))/(2.0*Sqrt(pi)))*((vesc**2.0/
&vmol_top(i)**2.0)+1.0)*Exp(-(vesc**2.0/vmol_top(i)**2.0))
      enddo

      !--------------------------------------------------------------------------
      ! ***  Jeans Escape plot ***
      !--------------------------------------------------------------------------

      !Writing out Jeans escape per molecule for plotting
      open(unit=75,file='jeans.dat',status='replace')
      do i=1,NMOLE
      	write(75,*) (trim(cmol(i))),Jeans(i)
      enddo
      close(75)


      end subroutine ESCAPE

