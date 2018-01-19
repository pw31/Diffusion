!-------------------------------------------------------------------------
      SUBROUTINE UPDATE_CRUST
!-------------------------------------------------------------------------
! ***  sets the "thickness" of the crust and its initial volume        ***
! ***  composition by considering solar abundances at the bottom of    ***
! ***  the atmosphere.                                                 ***
!-------------------------------------------------------------------------
      use GRID,ONLY: zz
      use STRUCT,ONLY: Temp,nHtot,nHeps,crust_depth,crust_beta,
     >                 crust_Neps,crust_Ncond,crust_gaseps
      use ELEMENTS,ONLY: NELEM,eps0,elnam
      use CHEMISTRY,ONLY: NMOLE,NELM,elnum,iel=>el
      use DUST_DATA,ONLY: NDUST,dust_vol,dust_nam,
     >                    dust_nel,dust_el,dust_nu
      use EXCHANGE,ONLY: inactive,nmol,chi
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real*8  :: Tg,nH,dz
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST),tmp(NELEM)
      real(kind=qp) :: sum_beta,Natmos,Ncrust,NtotH
      integer :: i,j,el,verbose

      Tg = Temp(1)
      nH = nHtot(1)
      dz = zz(2)-zz(1)
      eps0 = 0.Q0
      do i=1,NELM
        if (i==iel) cycle
        el = elnum(i)
        Natmos = nHeps(el,1)*dz
        Ncrust = crust_Neps(el)
        if (el==1) NtotH=Natmos+Ncrust
        eps0(el) = (Natmos+Ncrust)/NtotH
      enddo  
      inactive = .false.
      call EQUIL_COND(nH,Tg,eps,Sat,eldust,verbose)
      crust_gaseps = 0.Q0
      do i=1,NELM
        if (i==iel) cycle
        el = elnum(i)
        crust_gaseps(el) = eps(el) 
        nHeps(el,1) = nH*eps(el)
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
      !do i=1,NDUST
      !  if (eldust(i).le.0.Q0) cycle
      !  print*,dust_nam(i),REAL(crust_beta(i))
      !enddo  

      crust_Ncond(:) = 0.Q0                    ! crust column densities
      tmp(:) = 0.Q0
      do i=1,NDUST         
        if (eldust(i).le.0.Q0) cycle 
        crust_Ncond(i) = crust_beta(i)/dust_vol(i)
        do j=1,dust_nel(i)
          el = dust_el(i,j)
          tmp(el) = tmp(el) + dust_nu(i,j)*crust_Ncond(i)
        enddo
      enddo  
      do i=1,NELM
        if (i==iel) cycle
        el = elnum(i)
        if (tmp(el).le.0.Q0) cycle
        print*,elnam(el),REAL(crust_Neps(el)/tmp(el))
        crust_depth = crust_Neps(el)/tmp(el)   ! should all be the same
      enddo  
      crust_Ncond = crust_Ncond*crust_depth

      end
