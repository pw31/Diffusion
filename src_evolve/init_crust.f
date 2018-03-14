!-------------------------------------------------------------------------
      SUBROUTINE INIT_CRUST
!-------------------------------------------------------------------------
! ***  sets the "thickness" of the crust, its initial volume           ***
! ***  and column densities by considering phase equilibrium of a      ***
! ***  solar abundance gas at the bottom of the atmosphere.            ***
!-------------------------------------------------------------------------
      use PARAMETERS,ONLY: crust_kind,crust_thickness
      use GRID,ONLY: zz
      use STRUCT,ONLY: Temp,nHtot,press,crust_depth,crust_beta,
     >                 crust_Ncond,crust_Neps,crust_gaseps
      use ELEMENTS,ONLY: eps0,eps_solar,eps_meteor,eps_crust,elnam
      use CHEMISTRY,ONLY: NELEM,NELM,NMOLE,elnum,iel=>el
      use DUST_DATA,ONLY: NDUST,dust_vol,dust_nam,
     >                    dust_nel,dust_el,dust_nu
      use EXCHANGE,ONLY: inactive,nmol,chi,H,S,Fe,He
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real*8  :: Tg,nH,p,dz,PRESSURE
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST)
      real(kind=qp) :: sum_beta,fac
      integer :: i,j,el,it,verbose=0

      if (.not.allocated(nmol)) then
        allocate(nmol(NMOLE),chi(NDUST),inactive(NMOLE))
      endif  
      Tg = Temp(0)
      nH = nHtot(0)
      if (crust_kind==1) then
        eps0 = eps_solar
      else if (crust_kind==2) then 
        eps0 = eps_meteor
      else if (crust_kind==3) then 
        eps0 = eps_crust
      else if (crust_kind==4) then 
        !Tg = 500.d0
        eps0(:) = 1.E-99
        eps0(He) = 1.E+0 
        eps0(Fe) = 0.9E+0 
        eps0(S) = 1.E+0
      else
        stop "*** unknown crust_kind." 
      endif  
      inactive = .false.
      do it=1,99
        call EQUIL_COND(nH,Tg,eps,Sat,eldust,verbose)
        p = PRESSURE(Tg)
        print'(" iter=",I2," press=",2(1pE15.8)," eps0(H)=",1pE10.3)',
     >       it,p,press(0),eps0(H)
        fac = press(0)/p
        fac = MIN(3.0,MAX(0.3,fac))
        eps0 = eps0*fac
        if (ABS(1.d0-fac).lt.1.E-8) exit
      enddo  

      crust_gaseps = 0.Q0                      ! element abund. over crust
      print*
      print'(3x,2(A15))',"eps0","eps"
      do i=1,NELM
        if (i==iel) cycle
        el = elnum(i) 
        crust_gaseps(el) = eps(el) 
        print'(A3,2(1pE15.7))',elnam(el),eps0(el),eps(el)
      enddo 

      crust_beta(:) = 0.Q0                     ! crust volume composition
      sum_beta = 0.Q0
      do i=1,NDUST
        if (eldust(i).le.0.Q0) cycle
        crust_beta(i) = eldust(i)*dust_vol(i)
        sum_beta = sum_beta + crust_beta(i)
      enddo  
      crust_beta = crust_beta/sum_beta

      crust_Ncond(:) = 0.Q0                    ! condensed col.des. in crust
      crust_Neps(:)  = 0.Q0                    ! element col.dens. in crust
      crust_depth = crust_thickness
      do i=1,NDUST         
        if (eldust(i).le.0.Q0) cycle 
        crust_Ncond(i) = crust_depth*crust_beta(i)/dust_vol(i)
        do j=1,dust_nel(i)
          el = dust_el(i,j)
          crust_Neps(el) = crust_Neps(el) + dust_nu(i,j)*crust_Ncond(i)
        enddo
      enddo  

      print*
      print*," INIT_CRUST with crust_kind =",crust_kind
      print*,"            crust_depth[cm] =",crust_thickness

      end
