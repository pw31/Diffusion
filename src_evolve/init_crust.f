!-------------------------------------------------------------------------
      SUBROUTINE INIT_CRUST
!-------------------------------------------------------------------------
! ***  sets the "thickness" of the crust and its initial volume        ***
! ***  composition by considering solar abundances at the bottom of    ***
! ***  the atmosphere.                                                 ***
!-------------------------------------------------------------------------
      use NATURE,ONLY: km
      use GRID,ONLY: zz
      use STRUCT,ONLY: Temp,nHtot,crust_depth,
     >                 crust_beta,crust_eps,crust_gaseps
      use ELEMENTS,ONLY: NELEM,eps0,elnr,elnam,elcode
      use CHEMISTRY,ONLY: NMOLE
      use DUST_DATA,ONLY: NDUST,dust_vol,dust_nam,
     >                    dust_nel,dust_el,dust_nu
      use EXCHANGE,ONLY: inactive,nmol,chi
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real*8  :: Tg,nHges,delta_z
      real(kind=qp),dimension(NDUST) :: eps,Sat,eldust
      real(kind=qp) :: sum_beta,fac
      integer :: i,j,e,verbose

      allocate(nmol(NMOLE),chi(NDUST),inactive(NMOLE))
      Tg = Temp(1)
      nHges = nHtot(1)
      inactive = .false.
      call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
      crust_gaseps = 0.Q0
      do e=1,NELEM
        if (elcode(e)==0) cycle
        crust_gaseps(e) = eps(e) 
        print'(A3,2(1pE15.8))',elnam(e),eps0(e),eps(e)
      enddo  

      delta_z = zz(2)-zz(1)
      crust_depth   = 1.d0*km
      crust_eps(:)  = 0.Q0
      crust_beta(:) = 0.Q0
      sum_beta = 0.Q0
      do i=1,NDUST
        if (eldust(i).le.0.Q0) cycle
        crust_beta(i) = eldust(i)*dust_vol(i)
        sum_beta = sum_beta + crust_beta(i)
      enddo  
      crust_beta = crust_beta/sum_beta
      do i=1,NDUST
        if (eldust(i).le.0.Q0) cycle
        print*,dust_nam(i),REAL(crust_beta(i))
      enddo  

      do i=1,NDUST         ! crust abundances
        if (eldust(i).le.0.Q0) cycle 
        fac = crust_beta(i)/dust_vol(i)*(crust_depth/delta_z)/nHges
        do j=1,dust_nel(i)
          e = dust_el(i,j)
          crust_eps(e) = crust_eps(e) + dust_nu(i,j)*fac
        enddo
      enddo  
      eps0 = eps + crust_eps 
      call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
      do e=1,NELEM
        if (elcode(e)==0) cycle
        print'(A3,2(1pE15.8))',elnam(e),eps0(e),eps(e)
      enddo  

      end
