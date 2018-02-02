************************************************************************
      SUBROUTINE TEST
************************************************************************
      use NATURE,ONLY: bk,amu,mel,bar,km
      use ELEMENTS,ONLY: NELEM,NEPS,elcode,elnr,elnam,eps0,muH,mass
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_vol,dust_nel,dust_el,
     &                    dust_nu,qp
      use CHEMISTRY,ONLY: NMOLE,cmol,molmass
      use NUCLEATION,ONLY: NNUC,nuc,nuc_nam
      use REACTIONS,ONLY: NREAC,reac_sp
      use GRID,ONLY: Npoints
      use STRUCT,ONLY: nHtot,Temp,nHeps,rhoLj,rhoL3
      use EXCHANGE,ONLY: nel,nat,nmol,ipoint
      implicit none
      integer :: i,j,el,r,dustnr,stoich,NN
      real*8 :: nH,Tg,pp,rho,mu,sumn,sumnm,chinet,tt,dt
      real*8 :: Jst(NNUC),Nst(NNUC),chi(NDUST),bmix(NDUST)
      real*8 :: cr(NREAC),deps(NEPS)
      real*8,dimension(4+NDUST+NEPS) :: yy,FF
      real(kind=qp) :: eps(NELEM),Sat(NDUST)

      print*
      print*,"entering subroutine TEST ..."
      print*,"============================"

      do ipoint=50,50

        nH  = nHtot(ipoint)
        Tg  = Temp(ipoint)
        eps = eps0
        do i=1,NEPS
          el = elnr(i) 
          eps(el) = nHeps(i,ipoint)/nH
        enddo  

        !----------------------
        ! ***  call GGCHEM  ***
        !----------------------
        call GGCHEM(nH,Tg,eps,.false.,1)

        !----- compute pressure and mean molecular mass -----
        sumn  = nel
        sumnm = nel*mel
        do i=1,NELEM
          if (nat(i)==0.Q0) cycle
          sumn  = sumn  + nat(i)
          sumnm = sumnm + nat(i)*mass(i)
          !write(*,1000) trim(elnam(i)),nat(i),mass(i)/amu
        enddo  
        do i=1,NMOLE
          sumn  = sumn  + nmol(i)
          sumnm = sumnm + nmol(i)*molmass(i)
          !write(*,1000) trim(cmol(i)),nmol(i),molmass(i)/amu
        enddo
        mu = sumnm/sumn
        pp = sumn*bk*Tg
        print'(" n<H>[cm-3] =",1pE10.3,",  rho[g/cm3] =",1pE10.3)',
     >       nH,Tg
        print'("    mu[amu] =",0pF10.4,",      p[bar] =",1pE10.3)',
     >       mu/amu,pp/bar
        print*
        print*,'--- supersaturation ratios ---'
        call SUPERSAT(Tg,nat,nmol,Sat)
        do i=1,NDUST
          write(*,1100) trim(dust_nam(i)),Sat(i) 
        enddo  

        print*
        print*,'--- nucleation rates ---'
        call JSTAR(Tg,Sat,nat,nmol,Jst,Nst)
        do i=1,NNUC
          write(*,1200) trim(nuc_nam(i)),Sat(nuc(i)),Jst(i),Nst(i) 
        enddo  

        print*
        print*,'--- growth reactions ---'
        bmix = 1.0
        call GROWTH(Tg,Sat,bmix,nat,nmol,chinet,chi,cr)
        do i=1,NDUST
          write(*,2000) trim(dust_nam(i)),chi(i),chi(i)/dust_vol(i)
        enddo  

        print*
        print*,'--- element consumption ---'
        deps = 0.d0
        do r=1,NREAC
          dustnr = reac_sp(r,1,2)
          do j=1,dust_nel(dustnr)     ! loop over stoichiometry of condensate
            el = elcode(dust_el(dustnr,j))
            stoich = dust_nu(dustnr,j)
            deps(el) = deps(el) + stoich*cr(r)
          enddo  
        enddo
        do i=1,NEPS
          el = elnr(i) 
          write(*,3000),elnam(el),-deps(i)
        enddo  

        NN = 4+NDUST+NEPS
        yy(1:4) = rhoLj(0:3,ipoint)            ! rho*Lj    (j=0,1,2,3)
        yy(5:4+NDUST) = rhoL3(1:NDUST,ipoint)  ! rho*L3^s  (s=1...NDUST)
        do i=1,NEPS
          el = elnr(i) 
          yy(4+NDUST+i) = nHeps(i,ipoint)      ! element abundances
        enddo  
        tt = 0.d0
        call DUST_RHS(NN,tt,yy,FF,.true.) 
        dt = 9.e+99
        do j=4+NDUST+1,NN
          if (yy(j)>0.d0.and.FF(j).ne.0.d0) then
            dt = MIN(dt,0.01*yy(j)/ABS(FF(j)))
          endif  
        enddo  
        yy = yy + dt*FF
        call DUST_RHS(NN,tt,yy,FF,.true.)
        yy = yy + dt*FF
        print'(" yy=",99(1pE9.2))',yy
        print'(" FF=",99(1pE9.2))',FF
        print'(" dt=",1pE9.2)',dt
        call DUST_TSTEP(NN,yy,dt,1000.0)
        call DUST_RHS(NN,tt,yy,FF,.true.) 
      enddo
      stop

 1000 format(A12,1pE11.3,0pF9.3)
 1100 format(A16,"  S =",1pE10.3)
 1200 format(A10,"  S =",1pE10.3,"  Jst =",1pE10.3,
     &       " cm-3s-1  Nst =",1pE10.3)
 2000 format(A16,"  chi =",1pE11.3," cm/s,  ",1pE11.3," units/cm2/s")
 3000 format(A10,"  deps =",1pE11.3," atoms/cm2/s")
      end
