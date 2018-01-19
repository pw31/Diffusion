!-------------------------------------------------------------------------
      SUBROUTINE INIT_CRUST
!-------------------------------------------------------------------------
! ***  sets the "thickness" of the crust, its initial volume           ***
! ***  and column densities by considering phase equilibrium of a      ***
! ***  solar abundance gas at the bottom of the atmosphere.            ***
!-------------------------------------------------------------------------
      use NATURE,ONLY: km
      use GRID,ONLY: zz
      use STRUCT,ONLY: Temp,nHtot,crust_depth,crust_beta,crust_Ncond,
     >                 crust_Neps,crust_gaseps
      use ELEMENTS,ONLY: eps0,elnam
      use CHEMISTRY,ONLY: NELEM,NELM,NMOLE,elnum,iel=>el
      use DUST_DATA,ONLY: NDUST,dust_vol,dust_nam,
     >                    dust_nel,dust_el,dust_nu
      use EXCHANGE,ONLY: inactive,nmol,chi
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real*8  :: Tg,nH,dz
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST)
      real(kind=qp) :: sum_beta,fac
      integer :: i,j,el,verbose=0

      allocate(nmol(NMOLE),chi(NDUST),inactive(NMOLE))
      Tg = Temp(1)
      nH = nHtot(1)
      inactive = .false.
      call EQUIL_COND(nH,Tg,eps,Sat,eldust,verbose)

      crust_gaseps = 0.Q0                      ! element abund. over crust
      do i=1,NELM
        if (i==iel) cycle
        el = elnum(i) 
        crust_gaseps(el) = eps(el) 
        print'(A3,2(1pE15.8))',elnam(el),eps0(el),eps(el)
      enddo 

      crust_beta(:) = 0.Q0                     ! crust volume composition
      sum_beta = 0.Q0
      do i=1,NDUST
        if (eldust(i).le.0.Q0) cycle
        crust_beta(i) = eldust(i)*dust_vol(i)
        sum_beta = sum_beta + crust_beta(i)
      enddo  
      crust_beta = crust_beta/sum_beta

      crust_Ncond(:) = 0.Q0                       ! condensed col.des. in crust
      crust_Neps(:)  = 0.Q0                       ! element col.dens. in crust
      crust_depth = 1.d0*km                       ! initial thickness of crust
      crust_depth = 1.E-6
      crust_depth = 1.0
      do i=1,NDUST         
        if (eldust(i).le.0.Q0) cycle 
        crust_Ncond(i) = crust_depth*crust_beta(i)/dust_vol(i)
        do j=1,dust_nel(i)
          el = dust_el(i,j)
          crust_Neps(el) = crust_Neps(el) + dust_nu(i,j)*crust_Ncond(i)
        enddo
      enddo  

      end
