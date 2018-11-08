************************************************************************
      subroutine INIT_GRID
************************************************************************
      use GRID,ONLY: Npoints,zz
      use STRUCT,ONLY: Diff,nHtot
      use PARAMETERS,ONLY: Hp
      implicit none
      integer :: i,it
      integer :: Nsave=0,Nreach=3  !0,*=off, e.g. 3,3 to refine inner boundary
      real :: zmax

      !-----------------------------
      ! ***  set up test problem ***
      !-----------------------------
      allocate(zz(Npoints),Diff(Npoints),nHtot(Npoints))
      do i=1,Npoints-Nsave
        zz(i) = Hp*(REAL(i-1)/REAL(Npoints-Nsave-1))**1.0
      enddo

      !--------------------------------
      ! ***  refine inner boundary  ***
      !--------------------------------
      !print'(I2,99(1pE11.3))',0,zz(1:Nsave+2)
      do it=Nsave,1,-1
        zmax = zz(Nreach)
        do i=Npoints-it,Nreach,-1
          zz(i+1) = zz(i) 
        enddo
        do i=1,Nreach
          zz(i) = zmax*REAL(i-1)/REAL(Nreach)
        enddo
        !print'(I2,99(1pE11.3))',it,zz(1:Nsave+2)
      enddo
      !print'(I2,999(1pE11.3))',it,zz(1:Npoints)
      
      do i=1,Npoints
        Diff(i)  = 1.0                   ! diffusion coefficient
        nHtot(i) = 1.0                   ! density
        !nHtot(i) = 1.0*exp(-zz(i)/Hp)
      enddo  

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

 1000 format(A8,A11,A11,A11,A9,A12)
 1010 format(I4,I4,0pF11.3,1pE11.3,0pF11.2,0pF9.3,1pE12.3)
      end
