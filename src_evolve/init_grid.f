************************************************************************
      subroutine INIT_GRID
************************************************************************
      use NATURE,ONLY: bk,amu,km,bar
      use GRID,ONLY: N=>Npoints,zz,d1l,d1m,d1r,d2l,d2m,d2r
      use READMODEL,ONLY: Nlayers,zlay,Tlay,play,rholay,glay,Difflay
      use STRUCT,ONLY: Diff,rho,nHtot,Temp,press,mu,mols,atms,elec
      use CHEMISTRY,ONLY: NMOLE
      use ELEMENTS,ONLY: NELEM
      use PARAMETERS,ONLY: Hp,pmin,pmax
      implicit none
      integer :: i,j
      real*8 :: fac1,fac2,hmin,hmax,h1,h2,df
      real*8 :: hr(N),hl(N),f0(N),f1(N),f2(N)
      logical :: test=.false.

      allocate(zz(N),Diff(N),rho(N),nHtot(N),
     >         Temp(N),press(N),mu(N),
     >         mols(NMOLE,N),atms(NELEM,N),elec(N),
     >         d1l(1:N),d1m(1:N),d1r(1:N),
     >         d2l(1:N),d2m(1:N),d2r(1:N))

      hmin = 0.0
      hmax = 0.0
      do j=1,Nlayers
        !print*,j,zlay(j),play(j)/bar 
        if (play(j)>pmin*bar.and.hmax==0.0) hmax=zlay(j) 
        if (play(j)>pmax*bar.and.hmin==0.0) hmin=zlay(j) 
      enddo   
      !write(*,*) hmin,hmax
      do i=1,N
        zz(i) = hmin+(hmax-hmin)*(REAL(i-1)/REAL(N-1))**1.0
      enddo
  
      write(*,*)
      write(*,1000) "","z[km]","p[bar]","T[K]","mu[amu]","Diff[cm2/s]"
      j=Nlayers
      do i=1,N
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
        Temp(i)  = Tlay(j)*fac1 + Tlay(j-1)*fac2
        press(i) = EXP(LOG(play(j))*fac1 + LOG(play(j-1))*fac2)
        rho(i)   = EXP(LOG(rholay(j))*fac1 + LOG(rholay(j-1))*fac2)
        nHtot(i) = rho(i)/(1.4*amu)
        Diff(i)  = EXP(LOG(Difflay(j))*fac1 + LOG(Difflay(j-1))*fac2)
        mu(i)    = rho(i)/press(i)*bk*Temp(i)
        write(*,1010) i,j,(zz(i)-zz(1))/km,press(i)/bar,Temp(i),
     >                mu(i)/amu,Diff(i)
      enddo
      zz(:) = zz(:)-zz(1)

      !--- compute grid point differences ---
      do i=1,N-1
        hr(i) = zz(i+1)-zz(i)
      enddo  
      do i=2,N
        hl(i) = zz(i)-zz(i-1)
      enddo  

      !--- compute 1st and 2nd derivative coefficients ---
      do i=2,N-1
        d1l(i) = -hr(i)/((hr(i)+hl(i))*hl(i))
        d1m(i) =  (hr(i)-hl(i))/(hl(i)*hr(i))
        d1r(i) =  hl(i)/((hr(i)+hl(i))*hr(i)) 
        d2l(i) =  2.0/((hr(i)+hl(i))*hl(i))
        d2m(i) = -2.0/(hr(i)*hl(i))
        d2r(i) =  2.0/((hr(i)+hl(i))*hr(i))
      enddo

      !--- compute 1st derivative coefficients at boundaries ---
      !d1m(1) = -1.d0/hr(1)
      !d1r(1) =  1.d0/hr(1)  ! first order
      h1 = zz(2)-zz(1)
      h2 = zz(3)-zz(1)      
      d1l(1) = -(h1+h2)/(h1*h2)
      d1m(1) =  h2/(h1*(h2-h1))
      d1r(1) = -h1/(h2*(h2-h1))
      !d1l(N) = -1.d0/hl(N)
      !d1m(N) =  1.d0/hl(N)
      h1 = zz(N-1)-zz(N)
      h2 = zz(N-2)-zz(N)      
      d1r(N) = -(h1+h2)/(h1*h2)
      d1m(N) =  h2/(h1*(h2-h1))
      d1l(N) = -h1/(h2*(h2-h1))

      if (test) then
        !--- test derivatives ---
        do i=1,N
          f0(i) = 2.0*zz(i)**2-zz(i)+0.5
          f1(i) = 4.0*zz(i)-1.0
          f2(i) = 4.0
        enddo
        do i=2,N-1
          df = f0(i-1)*d1l(i) + f0(i)*d1m(i) + f0(i+1)*d1r(i)
          print'(I4,3(1pE14.6))',i,f1(i),df,f1(i)-df
        enddo
        do i=2,N-1
          df = f0(i-1)*d2l(i) + f0(i)*d2m(i) + f0(i+1)*d2r(i)
          print'(I4,3(1pE14.6))',i,f2(i),df,f2(i)-df
        enddo
        df = f0(1)*d1l(1) + f0(2)*d1m(1) + f0(3)*d1r(1) 
        print'(I4,3(1pE14.6))',1,f1(1),df,f1(1)-df
        df = f0(N-2)*d1l(N) + f0(N-1)*d1m(N) + f0(N)*d1r(N) 
        print'(I4,3(1pE14.6))',N,f1(N),df,f1(N)-df
        stop
      endif  

 1000 format(A8,A11,A11,A11,A9,A12)
 1010 format(I4,I4,0pF11.3,1pE11.3,0pF11.2,0pF9.3,1pE12.3)
      end
