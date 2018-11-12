************************************************************************
      subroutine INIT_GRID
************************************************************************
      use GRID,ONLY: N=>Npoints,zz
      use STRUCT,ONLY: Diff,nHtot
      use PARAMETERS,ONLY: Hp
      implicit none
      integer :: i,it
      integer :: Nsave=0,Nreach=3  !0,*=off, e.g. 3,3 to refine inner boundary
      real :: zmax

      !-----------------------------
      ! ***  set up test problem ***
      !-----------------------------
      allocate(zz(N),Diff(N),nHtot(N))
      do i=1,N-Nsave
        zz(i) = Hp*(REAL(i-1)/REAL(N-Nsave-1))**1.0
      enddo

      !--------------------------------
      ! ***  refine inner boundary  ***
      !--------------------------------
      !print'(I2,99(1pE11.3))',0,zz(1:Nsave+2)
      do it=Nsave,1,-1
        zmax = zz(Nreach)
        do i=N-it,Nreach,-1
          zz(i+1) = zz(i) 
        enddo
        do i=1,Nreach
          zz(i) = zmax*REAL(i-1)/REAL(Nreach)
        enddo
        !print'(I2,99(1pE11.3))',it,zz(1:Nsave+2)
      enddo
      !print'(I2,999(1pE11.3))',it,zz(1:N)

      !----------------------------------------------------
      ! ***  fill in density and diffusion coefficient  ***
      !----------------------------------------------------
      do i=1,N
        Diff(i)  = 1.0                   ! diffusion coefficient
        nHtot(i) = 1.0                   ! density
        !nHtot(i) = 1.0*exp(-zz(i)/Hp)
      enddo  

      call ZWEIGHTS

      end

************************************************************************
      subroutine ZWEIGHTS
************************************************************************
      use GRID,ONLY: N=>Npoints,zz,zweight
      implicit none
      integer :: i
      real,allocatable :: f(:)
      real :: h,dh,hm1,hm2,hp1,hp2,zmax,fint1,fint2
      logical :: test=.false.
      
      allocate(zweight(N),f(N))

      !-----------------------------------------------------------
      ! ***  3rd order z-integration weights, see Mathematica  ***
      !-----------------------------------------------------------
      zweight(:) = 0.0
      hm1 = zz(2)-zz(1)
      hp1 = zz(3)-zz(2)
      hp2 = zz(4)-zz(2)
      zweight(1) = zweight(1) + hm1*(3.0*hm1**2+6.0*hp1*hp2+4.0*hm1
     $     *(hp1+hp2))/(12.0*(hm1+hp1)*(hm1+hp2))
      zweight(2) = zweight(2) + hm1*(hm1**2+6.0*hp1*hp2+2.0*hm1*(hp1
     $     +hp2))/(12.0*hp1*hp2)
      zweight(3) = zweight(3) + hm1**3*(hm1+2.0*hp2)/(12.0*hp1*(hm1
     $     +hp1)*(hp1-hp2))
      zweight(4) = zweight(4) - hm1**3*(hm1+2.0*hp1)/(12.0*(hp1-hp2)
     $     *hp2*(hm1+hp2))
      hm1 = zz(N-2)-zz(N-3)
      hp1 = zz(N-1)-zz(N-2)
      hp2 = zz(N)  -zz(N-2)
      zweight(N-3) = zweight(N-3) - (hp1-hp2)**3*(hp1+hp2)/(12.0*hm1
     $     *(hm1+hp1)*(hm1+hp2)) 
      zweight(N-2) = zweight(N-2) + (hp1-hp2)**3*(2.0*hm1+hp1+hp2)/(12.0
     $     *hm1*hp1*hp2)
      zweight(N-1) = zweight(N-1) - (hp1-hp2)*(hp1*(4.0*hm1+3.0*hp1)+2.0
     $     *(hm1+hp1)*hp2+hp2**2)/(12.0*hp1*(hm1+hp1))
      zweight(N)   = zweight(N)   - (hp1-hp2)*(hp1**2+2.0*hp1*hp2+3.0
     $     *hp2**2+2.0*hm1*(hp1+2.0*hp2))/(12.0*hp2*(hm1+hp2))
      do i=2,N-2
        hm1 = zz(i)-zz(i-1)
        hp1 = zz(i+1)-zz(i)
        hp2 = zz(i+2)-zz(i)
        zweight(i-1) = zweight(i-1) + hp1**3*(hp1-2.0*hp2)/(12.0*hm1
     $       *(hm1+hp1)*(hm1+hp2))
        zweight(i)   = zweight(i)   -(hp1*(hp1*(2.0*hm1+hp1)-2.0*(3.0
     $       *hm1+hp1)*hp2))/(12.0*hm1*hp2)
        zweight(i+1) = zweight(i+1) +(hp1*(4.0*hm1*hp1+3.0*hp1**2-6.0
     $       *hm1*hp2-4.0*hp1*hp2))/(12.0*(hm1+hp1)*(hp1-hp2))
        zweight(i+2) = zweight(i+2) + hp1**3*(2.0*hm1+hp1)/(12.0*(hp1
     $       -hp2)*hp2*(hm1+hp2))
      enddo
      if (test) then
        h  = (zz(N)-zz(1))
        dh = h/(N-1)
        do i=1,N
          f(i) = 1.0 - zz(i)/h + (zz(i)/h)**2 - (zz(i)/h)**3
          print*,i,zweight(i)/dh
        enddo
        fint1 = (zz(N)-zz(1)) - 1.0/2.0/h*(zz(N)**2-zz(1)**2)
     &                     + 1.0/3.0/h**2*(zz(N)**3-zz(1)**3)
     &                     - 1.0/4.0/h**3*(zz(N)**4-zz(1)**4)
        fint2 = 0.0
        do i=1,N
          fint2 = fint2 + f(i)*zweight(i)
        enddo
        print*,fint1,fint2,fint2-fint1
        stop
      endif  
      end
      

************************************************************************
      subroutine FILLIN_GRID
************************************************************************
      use NATURE,ONLY: bk,amu,km,bar
      use GRID,ONLY: Npoints,zz
      use READMODEL,ONLY: Nlayers,zlay,Tlay,play,rholay,glay,Difflay
      use STRUCT,ONLY: Diff,rho,nHtot,T,press,mu
      use PARAMETERS,ONLY: Hp,pmin,pmax
      implicit none
      integer :: i,j
      real :: fac1,fac2,hmin,hmax

      allocate(zz(Npoints),Diff(Npoints),rho(Npoints),nHtot(Npoints),
     >         T(Npoints),press(Npoints),mu(Npoints))

      hmin = 0.0
      hmax = 0.0
      do j=1,Nlayers
        !print*,j,zlay(j),play(j)/bar 
        if (play(j)>pmin*bar.and.hmax==0.0) hmax=zlay(j) 
        if (play(j)>pmax*bar.and.hmin==0.0) hmin=zlay(j) 
      enddo   
      write(*,*) Hp,hmin,hmax
      do i=1,Npoints
        zz(i) = hmin+(hmax-hmin)*(REAL(i-1)/REAL(Npoints-1))**1.0
      enddo
  
      write(*,*)
      write(*,1000) "","z[km]","p[bar]","T[K]","mu[amu]","Diff[cm2/s]"
      j=Nlayers
      do i=1,Npoints
        do 
          if (zz(i)<zlay(j-1).or.j==2) exit
          j=j-1
        enddo  
        !--- check zlay(j)<=zz(i)<zlay(j-1) ---
        if (zlay(j)>zz(i).or.zz(i)>zlay(j-1)) then
          print*,"could not bracket z-gridpoint" 
          print*,i,j,zlay(j),zz(i),zlay(j-1)
          stop
        endif   
        fac2 = (zlay(j)-zz(i))/(zlay(j)-zlay(j-1))
        fac1 = 1.0-fac2
        !print*,i,j,zlay(j),zz(i),zlay(j-1),fac2
        T(i)     = Tlay(j)*fac1 + Tlay(j-1)*fac2
        press(i) = EXP(LOG(play(j))*fac1 + LOG(play(j-1))*fac2)
        rho(i)   = EXP(LOG(rholay(j))*fac1 + LOG(rholay(j-1))*fac2)
        nHtot(i) = rho(i)/(1.4*amu)
        Diff(i)  = EXP(LOG(Difflay(j))*fac1 + LOG(Difflay(j-1))*fac2)
        mu(i)    = rho(i)/press(i)*bk*T(i)
        write(*,1010) i,j,(zz(i)-zz(1))/km,press(i)/bar,T(i),
     >                mu(i)/amu,Diff(i)
      enddo
      zz(:) = zz(:)-zz(1)

      call ZWEIGHTS
      
 1000 format(A8,A11,A11,A11,A9,A12)
 1010 format(I4,I4,0pF11.3,1pE11.3,0pF11.2,0pF9.3,1pE12.3)
      end
