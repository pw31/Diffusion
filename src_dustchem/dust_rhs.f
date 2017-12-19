************************************************************************
      SUBROUTINE DUST_RHS(NN,t,yy,FF,outp,ifail)
************************************************************************
      use PARAMETERS,ONLY: Vl
      use NATURE,ONLY: bk,amu,mel,bar,km,pi,mic
      use ELEMENTS,ONLY: NELEM,NEPS,elcode,elnr,elnam,eps0,muH,mass
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_vol,dust_nel,dust_el,
     &                    dust_nu,dust_vol,dust_rho,qp
      use CHEMISTRY,ONLY: NMOLE,cmol,molmass
      use NUCLEATION,ONLY: NNUC,nuc,nuc_nam
      use REACTIONS,ONLY: NREAC,reac_sp
      use STRUCT,ONLY: nHtot,Temp,press
      use EXCHANGE,ONLY: nel,nat,nmol,ipoint,Jst,chi
      implicit none
      integer,intent(IN) :: NN
      real*8,intent(IN)  :: t,yy(NN)
      real*8,intent(OUT) :: FF(NN)
      logical,intent(in) :: outp
      integer,intent(out) :: ifail
      integer :: i,j,s,el,r,dustnr,stoich,ierr
      real*8 :: nH,Tg,pp,rhod,bsum,Jsum,chinet,V0,Nl,amean,rhoL(0:3)
      real*8 :: Nst(NNUC),bmix(NDUST),cr(NREAC)
      real(kind=qp) :: eps(NELEM),Sat(NDUST)
      character(len=9999) :: info
      logical :: IS_NAN

      !-----------------------------------------------
      ! ***  get atmospheric properties at ipoint  ***
      !-----------------------------------------------
      Tg  = Temp(ipoint)
      nH  = nHtot(ipoint)
      pp  = press(ipoint)
      ifail = 0

      !--------------------------------------------------------------
      ! ***  extract element abundances and dust moments from yy  ***
      !--------------------------------------------------------------
      rhoL(0:3) = yy(1:4)                 ! dust moments rho*Lj (j=0...3)
      bmix(1:NDUST) = yy(5:4+NDUST)       ! dust volume fractions
      eps = eps0
      do i=1,NEPS
        el = elnr(i) 
        eps(el) = yy(4+NDUST+i)/nH        ! element abundances
        if ((eps(el).lt.0.Q0).or.IS_NAN(REAL(eps(el)))) then
          print*,"*** eps=NaN in dust_rhs"
          print*,ipoint,nH,Tg
          print*,yy
          ifail = -1
          stop
          return
        endif
      enddo

      !----------------------
      ! ***  call GGCHEM  ***
      !----------------------
      call GGCHEM(nH,Tg,eps,.false.,0)

      !------------------------------------------------------
      ! ***  supersaturation ratios and nucleation rates  ***
      !------------------------------------------------------
      call SUPERSAT(Tg,nat,nmol,Sat)
      call JSTAR(Tg,Sat,nat,nmol,Jst,Nst)
      Jsum = 0.d0
      do i=1,NNUC
        Jsum = Jsum + Jst(i)
      enddo 

      !------------------------------------------
      ! ***  extract dust quantities from yy  ***
      !------------------------------------------
      bsum = 0.d0                         
      rhod = 0.d0                         ! dust material density
      do i=1,NDUST
        bsum = bsum + bmix(i)
        rhod = rhod + bmix(i)*dust_rho(i)
      enddo
      if (bsum>0.d0) then
        !print*,"check rhoL3:",bsum,rhoL(3)
        bmix = bmix/bsum
        rhod = rhod/bsum
      else
        bsum = 0.d0                         
        bmix = 0.d0 
        rhod = 0.d0  
        do i=1,NNUC
          s = nuc(i)
          V0 = dust_vol(s)
          bmix(s) = Jst(i)*V0
          bsum = bsum + bmix(s)
          rhod = rhod + bmix(s)*dust_rho(s)
        enddo
        if (bsum>0.d0) then
          bmix = bmix/bsum
          rhod = rhod/bsum
        else
          rhod = 3.d0 
        endif   
      endif  

      !---------------------------
      ! ***  growth reactions  ***
      !---------------------------
      call GROWTH(Tg,Sat,bmix,nat,nmol,chinet,chi,cr)

      !-----------------------------
      ! ***  compute rhs vector  ***
      !-----------------------------
      FF(:) = 0.d0
      FF(1) = Jsum
      do i=1,3
        FF(i+1) = Vl**(DBLE(i)/3.d0)*Jsum                        ! nucleation
     &          + DBLE(i)/3.d0*chinet*rhoL(i-1)                  ! net growth
      enddo  
      do i=1,NNUC
        s  = nuc(i)
        V0 = dust_vol(s)
        Nl = Vl/V0
        FF(4+s) = FF(4+s) + Vl*Jst(i)                            ! nucleation
        do j=1,dust_nel(s)            
          el = elcode(dust_el(s,j))
          stoich = dust_nu(s,j)
          FF(4+NDUST+el) = FF(4+NDUST+el) - stoich*Nl*Jst(i)     ! consumption
        enddo   
      enddo  
      do s=1,NDUST
        FF(4+s) = FF(4+s) + chi(s)*rhoL(2)                       ! net growth
      enddo  
      do r=1,NREAC
        s = reac_sp(r,1,2)
        do j=1,dust_nel(s)                        
          el = elcode(dust_el(s,j))
          stoich = dust_nu(s,j)
          FF(4+NDUST+el) = FF(4+NDUST+el) - stoich*cr(r)*rhoL(2) ! consumption
        enddo  
      enddo  

      if (outp) then
        print* 
        print'(I4,"  n<H>=",1pE9.2,"  p/bar=",1pE9.2,"  T/K=",0pF8.2)',
     >        ipoint,nH,pp/bar,Tg 
        do i=1,NEPS
          el = elnr(i) 
          print'(A4,2(1pE11.3))',elnam(el),eps(el),eps(el)/eps0(el) 
        enddo  
        if (rhoL(0)>0.d0) then
          amean = (3.d0/(4.d0*pi))**(1.d0/3.d0) * rhoL(1)/rhoL(0)
        else
          amean = 0.d0 
        endif
        print'("  rho*L=",4(1pE10.2),"  <a>/mic=",1pE10.2)',
     >        rhoL(0:3),amean/mic
        print'("    Jst=",99(1pE10.2))',Jst
        print'(8x,10x,6x,99(A10))',(trim(dust_nam(i)),i=1,NDUST) 
        print'(        8x,    10x," bmix=",99(1pE10.2))',bmix
        print'(        8x,    10x,"  Sat=",99(1pE10.2))',Sat
        print'(" chinet=",1pE10.2,"  chi=",99(1pE10.2))',chinet,chi
      endif  

#ifdef CHECK_NAN
      do i=1,NN
        if (IS_NAN(FF(i))) then
          print*,"*** DUST_RHS returns NaN",ipoint
          print*,"Jst=",Jst
          print*,"chi=",chi
          print*,"FF=",FF
          call TRACEBACKQQ(info,ierr)
          print*,ierr,info
          stop "*** stop"
        endif  
      enddo
#endif

      end



