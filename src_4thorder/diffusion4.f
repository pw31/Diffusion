************************************************************************
      subroutine DIFFUSION_EXPLICIT
************************************************************************
      use NATURE,ONLY: pi
      use GRID,ONLY: N=>Npoints,zz,xx,d1l2,d1l1,d1m,d1r1,
     >d1r2,d2l2,d2l1,d2m,d2r1,d2r2
      use PARAMETERS,ONLY: Hp,tnull,init,bc_low,bc_high,
     >                     influx,outflux,inrate,outrate,vin,vout,
     >                     Nout,outtime
      use STRUCT,ONLY: Diff,nHtot
      implicit none
      integer :: i,it
      real*8,dimension(N) :: rate
      real*8 :: D,nD,d1,d2,d1nD,jdiff,time,dt,tend
      real*8 :: ww,AA,z0,x1ana,xNana
      real*8 :: ntot,ntot2
      character :: CR = CHAR(13)

      !-------------------
      ! ***  timestep  ***
      !-------------------
      dt = 9.D+99
      do i=2,N
        D  = 0.5*(Diff(i-1)+Diff(i)) 
        dt = MIN(dt,0.33*(zz(i)-zz(i-1))**2/D)
      enddo  
      dt = dt*(1.d0+1.d-12)
      print'("explicit timestep =",1pE10.2," s")',dt
      print*
      
      ntot = 0.d0
      do i=1,N-1
        ntot = ntot + 0.5*(nHtot(i)*xx(i)+nHtot(i+1)*xx(i+1))
     >                   *(zz(i+1)-zz(i))
      enddo  
      
      time = tnull
      tend = outtime(Nout)
      Nout = 1
      open(unit=1,file='out.dat',status='replace')
      write(1,*) N,init 
      write(1,'(9999(1pE16.8))') zz
      write(1,'("time[s]=",1pE12.5)') time 
      write(1,'(9999(1pE16.8))') (MAX(xx(i),1.E-99),i=1,N)

      do it=1,9999999

        !------------------------------
        ! ***  boundary conditions  ***
        !------------------------------
        if (init==4) then  ! the analytic solution has flux=flux(t)
          z0 = 0.5*(zz(N)+zz(1))
          !ww = 2.d0*SQRT(nHtot(1)*Diff(1)*time)
          !AA = SQRT(tnull/time)
          !x1ana = AA*exp(-((zz(1)-z0)/ww)**2)
          !xNana = AA*exp(-((zz(N)-z0)/ww)**2)
          !influx  = 0.5*(zz(1)-z0)*x1ana/time   ! analytic solution
          !outflux = 0.5*(zz(N)-z0)*xNana/time   ! for fluxes
          influx  = 0.5*(zz(1)-z0)*xx(1)/time
          outflux = 0.5*(zz(N)-z0)*xx(N)/time
        endif  
        nD = nHtot(1)*Diff(1)  
        if (bc_low==1) then                      ! fixed concentration
          influx = -nD*( d1l2(1)*xx(1)
     >                  +d1l1(1)*xx(2)
     >                  + d1m(1)*xx(3) 
     >                  +d1r1(1)*xx(4)
     >                  +d1r2(1)*xx(5))
        else if (bc_low==2) then                 ! fixed flux
          xx(1) = (-influx/nD -d1l1(1)*xx(2)
     >                        -d1m(1) *xx(3)
     >                        -d1r1(1)*xx(4)
     >                        -d1r2(1)*xx(5) )/d1l2(1)
        else if (bc_low==3) then   
          influx = nHtot(1)*xx(1)*inrate*vin
          xx(1) = -( d1l1(1)*xx(2)               ! fixed inrate
     >              +d1m(1) *xx(3)
     >              +d1r1(1)*xx(4)
     >              +d1r2(1)*xx(5) )/(d1l2(1) + inrate*vin/Diff(1))
        endif  
        nD = nHtot(N)*Diff(N)  
        if (bc_high==1) then                     ! fixed concentration
          outflux = -nD*( d1l2(N)*xx(N-4)
     >                   +d1l1(N)*xx(N-3)
     >                   +d1m(N) *xx(N-2)
     >                   +d1r1(N)*xx(N-1)
     >                   +d1r2(N)*xx(N)  )
        else if (bc_high==2) then                ! fixed outrate
          xx(N) = (-outflux/nD -d1l2(N)*xx(N-4)
     >                         -d1l1(N)*xx(N-3) 
     >                         -d1m(N) *xx(N-2)
     >                         -d1r1(N)*xx(N-1) )/d1r2(N)
        else if (bc_high==3) then   
          outflux = nHtot(N)*xx(N)*outrate*vout  
          xx(N) = -( d1l2(N)*xx(N-4)
     >              +d1l1(N)*xx(N-3)
     >              +d1m(N) *xx(N-2) 
     >              +d1r1(N)*xx(N-1) )/(d1r2(N) + outrate*vout/Diff(N))
        endif   
        ntot = ntot + (influx - outflux)*0.5*dt
        
        !-------------------------------------------
        ! ***  d/dt(nH*x) = d/dz(nH*Diff*dx/dz)  ***
        !-------------------------------------------
        rate(:) = 0.0
        do i=3,N-2
          nD   = nHtot(i)*Diff(i)  
          d1   = d1l2(i)*xx(i-2)
     >         + d1l1(i)*xx(i-1)
     >         + d1m(i) *xx(i) 
     >         + d1r1(i)*xx(i+1)
     >         + d1r2(i)*xx(i+2)
          d2   = d2l2(i)*xx(i-2)
     >         + d2l1(i)*xx(i-1)
     >         + d2m(i) *xx(i) 
     >         + d2r1(i)*xx(i+1)
     >         + d2r2(i)*xx(i+2)
          d1nD = d1l2(i)*nHtot(i-2)*Diff(i-2)
     >         + d1l1(i)*nHtot(i-1)*Diff(i-1)
     >         + d1m(i) *nHtot(i)  *Diff(i) 
     >         + d1r1(i)*nHtot(i+1)*Diff(i+1)
     >         + d1r2(i)*nHtot(i+2)*Diff(i+2)
          rate(i) = nD*d2 + d1nD*d1
        enddo
        !--- i=2 ---
        nD   = nHtot(2)*Diff(2)  
        d1   = d1l2(2)*xx(1)
     >       + d1l1(2)*xx(2)
     >       + d1m(2) *xx(3) 
     >       + d1r1(2)*xx(4)
     >       + d1r2(2)*xx(5)
        d2   = d2l2(2)*xx(1)
     >       + d2l1(2)*xx(2)
     >       + d2m(2) *xx(3) 
     >       + d2r1(2)*xx(4)
     >       + d2r2(2)*xx(5)
        d1nD = d1l2(2)*nHtot(1)*Diff(1)
     >       + d1l1(2)*nHtot(2)*Diff(2)
     >       + d1m(2) *nHtot(3)*Diff(3) 
     >       + d1r1(2)*nHtot(4)*Diff(4)
     >       + d1r2(2)*nHtot(5)*Diff(5)
        rate(2) = nD*d2 + d1nD*d1
        !--- i=N-1 ---
        nD   = nHtot(N-1)*Diff(N-1)  
        d1   = d1l2(N-1)*xx(N-4)
     >       + d1l1(N-1)*xx(N-3)
     >       + d1m(N-1) *xx(N-2) 
     >       + d1r1(N-1)*xx(N-1)
     >       + d1r2(N-1)*xx(N)
        d2   = d2l2(N-1)*xx(N-4)
     >       + d2l1(N-1)*xx(N-3)
     >       + d2m(N-1) *xx(N-2) 
     >       + d2r1(N-1)*xx(N-1)
     >       + d2r2(N-1)*xx(N)
        d1nD = d1l2(N-1)*nHtot(N-4)*Diff(N-4)
     >       + d1l1(N-1)*nHtot(N-3)*Diff(N-3)
     >       + d1m(N-1) *nHtot(N-2)*Diff(N-2) 
     >       + d1r1(N-1)*nHtot(N-1)*Diff(N-1)
     >       + d1r2(N-1)*nHtot(N)  *Diff(N)
        rate(N-1) = nD*d2 + d1nD*d1

        !----------------------------
        ! ***  explicit timestep  ***
        !----------------------------
        do i=2,N-1
          xx(i) = xx(i) + rate(i)*dt/nHtot(i)
        enddo  
        time = time + dt
        write(*,'(TL10,I8,A,$)') it,CR

        !------------------------------
        ! ***  boundary conditions  ***
        !------------------------------
        if (init==4) then  ! the analytic solution has flux=flux(t)
          z0 = 0.5*(zz(N)+zz(1))
          !ww = 2.d0*SQRT(nHtot(1)*Diff(1)*time)
          !AA = SQRT(tnull/time)
          !x1ana = AA*exp(-((zz(1)-z0)/ww)**2)
          !xNana = AA*exp(-((zz(N)-z0)/ww)**2)
          !influx  = 0.5*(zz(1)-z0)*x1ana/time   ! analytic solution
          !outflux = 0.5*(zz(N)-z0)*xNana/time   ! for fluxes
          influx  = 0.5*(zz(1)-z0)*xx(1)/time
          outflux = 0.5*(zz(N)-z0)*xx(N)/time
        endif  
        nD = nHtot(1)*Diff(1)  
        if (bc_low==1) then                      ! fixed concentration
          influx = -nD*( d1l2(1)*xx(1)
     >                  +d1l1(1)*xx(2)
     >                  + d1m(1)*xx(3) 
     >                  +d1r1(1)*xx(4)
     >                  +d1r2(1)*xx(5))
        else if (bc_low==2) then                 ! fixed flux
          xx(1) = (-influx/nD -d1l1(1)*xx(2)
     >                        -d1m(1) *xx(3)
     >                        -d1r1(1)*xx(4)
     >                        -d1r2(1)*xx(5) )/d1l2(1)
        else if (bc_low==3) then   
          influx = nHtot(1)*xx(1)*inrate*vin
          xx(1) = -( d1l1(1)*xx(2)               ! fixed inrate
     >              +d1m(1) *xx(3)
     >              +d1r1(1)*xx(4)
     >              +d1r2(1)*xx(5) )/(d1l2(1) + inrate*vin/Diff(1))
        endif  
        nD = nHtot(N)*Diff(N)  
        if (bc_high==1) then                     ! fixed concentration
          outflux = -nD*( d1l2(N)*xx(N-4)
     >                   +d1l1(N)*xx(N-3)
     >                   +d1m(N) *xx(N-2)
     >                   +d1r1(N)*xx(N-1)
     >                   +d1r2(N)*xx(N)  )
        else if (bc_high==2) then                ! fixed outrate
          xx(N) = (-outflux/nD -d1l2(N)*xx(N-4)
     >                         -d1l1(N)*xx(N-3) 
     >                         -d1m(N) *xx(N-2)
     >                         -d1r1(N)*xx(N-1) )/d1r2(N)
        else if (bc_high==3) then   
          outflux = nHtot(N)*xx(N)*outrate*vout  
          xx(N) = -( d1l2(N)*xx(N-4)
     >              +d1l1(N)*xx(N-3)
     >              +d1m(N) *xx(N-2) 
     >              +d1r1(N)*xx(N-1) )/(d1r2(N) + outrate*vout/Diff(N))
        endif   
        ntot = ntot + (influx - outflux)*0.5*dt
        
        if (time>outtime(Nout)) then
          write(*,'(I8," output t=",1pE14.6," s")') it,time
          write(1,'("time[s]=",1pE12.5)') time 
          write(1,'(9999(1pE16.8))') (MAX(xx(i),1.E-99),i=1,N)
          Nout = Nout + 1
          ntot2 = 0.0
          do i=1,N-1
            ntot2 = ntot2 + 0.5*(nHtot(i)*xx(i)+nHtot(i+1)*xx(i+1))
     >                         *(zz(i+1)-zz(i))
          enddo  
          print'("  total=",2(1pE14.6)," , dev=",0pF8.5,"%")',
     >         ntot,ntot2,(ntot/ntot2-1.0)*100.0
        endif  

        if (time>tend) exit

      enddo  

      print*
      print'(A4,99(A12))','grid','zz/Hp','eps','jdiff'
      do i=1,N
        nD = nHtot(i)*Diff(i)
        if (i==1) then
          d1 = d1l2(1)*xx(1)
     >        +d1l1(1)*xx(2)
     >        +d1m(1) *xx(3)
     >        +d1r1(1)*xx(4)
     >        +d1r2(1)*xx(5)
        else if (i==2) then
          d1 = d1l2(2)*xx(1)
     >        +d1l1(2)*xx(2)
     >        +d1m(2) *xx(3)
     >        +d1r1(2)*xx(4)
     >        +d1r2(2)*xx(5)
        else if (i==N-1) then
          d1 = d1l2(N-1)*xx(N-4)
     >        +d1l1(N-1)*xx(N-3)
     >        +d1m(N-1) *xx(N-2)
     >        +d1r1(N-1)*xx(N-1)
     >        +d1r2(N-1)*xx(N)
        else if (i==N) then
          d1 = d1l2(N)*xx(N-4)
     >        +d1l1(N)*xx(N-3)
     >        +d1m(N) *xx(N-2)
     >        +d1r1(N)*xx(N-1)
     >        +d1r2(N)*xx(N)
        else
          d1 = d1l2(i)*xx(i-2)
     >        +d1l1(i)*xx(i-1)
     >        +d1m(i) *xx(i)
     >        +d1r1(i)*xx(i+1)
     >        +d1r2(i)*xx(i+2)
        endif
        jdiff = -nD*d1
        print'(I4,0pF12.5,99(1pE12.4))',i,zz(i)/Hp,xx(i),jdiff
      enddo

      end

************************************************************************
      subroutine DIFFUSION_IMPLICIT
************************************************************************
c      use NATURE,ONLY: pi
c      use GRID,ONLY: N=>Npoints,zz,xx,d1l,d1m,d1r,d2l,d2m,d2r
c      use PARAMETERS,ONLY: Hp,bc_low,bc_high,init,tnull,
c     >                     influx,outflux,inrate,outrate,vin,vout,
c     >                     Nout,outtime,tfac
c      use STRUCT,ONLY: Diff,nHtot
c      implicit none
c      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
c      integer :: i,j,k,it,ipvt(N),info
c      real*8,dimension(N,N) :: A,B,sum
c      real(kind=qp),dimension(N,N) :: Awork
c      real(kind=qp) :: det(2),work(N)
c      real*8,dimension(N) :: xnew,rest
c      real*8 :: D,nD,d1,d1nD,jdiff,time,dt,tend
c      real*8 :: z0,ww,AA,x1ana,xNana,ntot,ntot2,err
c      logical :: check=.true.
c      character :: CR = CHAR(13)
c
c      !-------------------
c      ! ***  timestep  ***
c      !-------------------
c      dt = 9.D+99
c      do i=2,N
c        D  = 0.5*(Diff(i-1)+Diff(i))
c        dt = MIN(dt,0.5*(zz(i)-zz(i-1))**2/D)
c      enddo
c      print*,"explicit timestep=",dt
c      dt = dt*tfac
c      dt = dt*(1.d0+1.d-12)
c      print*," applied timestep=",dt
c
c      !-----------------------------
c      ! ***  fill in big matrix  ***
c      !-----------------------------
c      A(:,:) = 0.d0
c      do i=2,N-1
c        nD   = nHtot(i)*Diff(i)
c        d1nD = d1l(i)*nHtot(i-1)*Diff(i-1)
c     >       + d1m(i)*nHtot(i)  *Diff(i)
c     >       + d1r(i)*nHtot(i+1)*Diff(i+1)
c        A(i,i-1) = A(i,i-1) - dt*( d1nD*d1l(i) + nD*d2l(i) )
c        A(i,i)   = A(i,i)   - dt*( d1nD*d1m(i) + nD*d2m(i) )
c        A(i,i+1) = A(i,i+1) - dt*( d1nD*d1r(i) + nD*d2r(i) )
c      enddo
c      do i=1,N
c        A(i,:) = A(i,:)/nHtot(i)    ! unitless
c      enddo
c      !--------------------------
c      ! ***  add unit matrix  ***
c      !--------------------------
c      do i=1,N
c        A(i,i) = A(i,i) + 1.d0
c      enddo
c      !------------------------------
c      ! ***  boundary conditions  ***
c      !------------------------------
c      if (bc_low==1) then
c        !--- nothing to do
c      else if (bc_low==2) then
c        A(1,1) = 1.d0
c        A(1,2) = d1m(1)/d1l(1)
c        A(1,3) = d1r(1)/d1l(1)
c      else if (bc_low==3) then
c        A(1,1) = 1.d0+inrate*vin/Diff(1)/d1l(1)
c        A(1,2) = d1m(1)/d1l(1)
c        A(1,3) = d1r(1)/d1l(1)
c      endif
c      if (bc_high==1) then
c        !--- nothing to do
c      else if (bc_high==2) then
c        A(N,N-2) = d1l(N)/d1r(N)
c        A(N,N-1) = d1m(N)/d1r(N)
c        A(N,N)   = 1.d0
c      else if (bc_high==3) then
c        A(N,N-2) = d1l(N)/d1r(N)
c        A(N,N-1) = d1m(N)/d1r(N)
c        A(N,N)   = 1.d0+outrate*vout/Diff(N)/d1r(N)
c      endif
c
c      !------------------------
c      ! ***  invert matrix  ***
c      !------------------------
c      Awork = A
c      call QGEFA ( Awork, N, N, ipvt, info )
c      call QGEDI ( Awork, N, N, ipvt, det, work, 1 )
c      B = Awork
c      if (info.ne.0) then
c        print*,"*** singular matrix in QGEFA: info=",info
c        stop
c      endif
c      if (check) then
c        do i=1,N
c          write(99,'(9999(1pE11.3))') (A(i,j),j=1,N)
c        enddo
c        write(99,*)
c        do i=1,N
c          write(99,'(9999(1pE11.3))') (B(i,j),j=1,N)
c        enddo
c        !--- test A*B=1 ---
c        err = 0.d0
c        write(99,*)
c        do i=1,N
c          do j=1,N
c            sum(i,j) = 0.d0
c            do k=1,N
c              sum(i,j) = sum(i,j) + A(i,k)*B(k,j)
c            enddo
c            if (i==j) then
c              err = max(err,ABS(sum(i,j)-1.d0))
c            else
c              err = max(err,ABS(sum(i,j)))
c            endif
c          enddo
c          write(99,'(9999(1pE11.3))') (sum(i,j),j=1,N)
c        enddo
c        write(99,*) "maximum error=",err
c        print*,"matrix inversion error=",err
c      endif
c
c      ntot = 0.0
c      do i=1,N-1
c        ntot = ntot + 0.5*(nHtot(i)*xx(i)+nHtot(i+1)*xx(i+1))
c     >                   *(zz(i+1)-zz(i))
c      enddo
c
c      time = tnull
c      tend = outtime(Nout)
c      Nout  = 1
c      print*
c      open(unit=1,file='out.dat',status='replace')
c      write(1,*) N,init
c      write(1,'(9999(1pE16.8))') zz
c      write(1,'("time[s]=",1pE12.5)') time
c      write(1,'(9999(1pE16.8))') (MAX(xx(i),1.E-99),i=1,N)
c
c      do it=1,9999999
c
c        !------------------------------
c        ! ***  boundary conditions  ***
c        !------------------------------
c        if (init==4) then  ! the analytic solution has flux=flux(t)
c          z0 = 0.5*(zz(N)+zz(1))
c          !ww = 2.d0*SQRT(nHtot(1)*Diff(1)*time)
c          !AA = SQRT(tnull/time)
c          !x1ana = AA*exp(-((zz(1)-z0)/ww)**2)
c          !xNana = AA*exp(-((zz(N)-z0)/ww)**2)
c          !influx  = 0.5*(zz(1)-z0)*x1ana/time   ! analytic solution
c          !outflux = 0.5*(zz(N)-z0)*xNana/time   ! for fluxes
c          influx  = 0.5*(zz(1)-z0)*xx(1)/time
c          outflux = 0.5*(zz(N)-z0)*xx(N)/time
c        endif
c        rest = xx
c        nD = nHtot(1)*Diff(1)
c        if (bc_low==1) then
c          influx  = -nD*(d1l(1)*xx(1) + d1m(1)*xx(2) + d1r(1)*xx(3))
c        else if (bc_low==2) then
c          rest(1) = -influx/nD/d1l(1)
c        else if (bc_low==3) then
c          influx  = inrate*nHtot(1)*xx(1)*vin
c          rest(1) = 0.d0
c        endif
c        nD = nHtot(N)*Diff(N)
c        if (bc_high==1) then
c          outflux = -nD*(d1l(N)*xx(N-2) + d1m(N)*xx(N-1) + d1r(N)*xx(N))
c        else if (bc_high==2) then
c          rest(N) = -outflux/nD/d1r(N)
c        else if (bc_high==3) then
c          outflux = outrate*nHtot(N)*xx(N)*vout
c          rest(N) = 0.d0
c        endif
c        !ntot = ntot + (influx-outflux)*dt
c
c        !-----------------------
c        ! ***  source terms  ***
c        !-----------------------
c        !rest(N/3) = rest(N/3) + 500.0*dt
c
c        !----------------------------
c        ! ***  implicit timestep  ***
c        !----------------------------
c        do i=1,N
c          xnew(i) = 0.d0
c          do j=1,N
c            xnew(i) = xnew(i) + B(i,j)*rest(j)
c          enddo
c        enddo
c        xx(:) = xnew(:)
c        nD = nHtot(1)*Diff(1)
c        if (bc_low==1) then
c          influx  = -nD*(d1l(1)*xx(1) + d1m(1)*xx(2) + d1r(1)*xx(3))
c        else if (bc_low==3) then
c          influx  = inrate*nHtot(1)*xx(1)*vin
c        endif
c        nD = nHtot(N)*Diff(N)
c        if (bc_high==1) then
c          outflux = -nD*(d1l(N)*xx(N-2) + d1m(N)*xx(N-1) + d1r(N)*xx(N))
c        else if (bc_high==3) then
c          outflux = outrate*nHtot(N)*xx(N)*vout
c        endif
c        ntot = ntot + (influx-outflux)*dt
c        time = time + dt
c        write(*,'(TL10,I8,A,$)') it,CR
c
c        if (time>outtime(Nout)) then
c          write(*,'(I8," output t=",1pE11.3," s")') it,time
c          write(1,'("time[s]=",1pE12.5)') time
c          write(1,'(9999(1pE16.8))') (MAX(xx(i),1.E-99),i=1,N)
c          Nout = Nout + 1
c        endif
c
c        if (time>tend) exit
c
c      enddo
c
c      print*
c      print'(A4,99(A12))','grid','zz/Hp','eps','jdiff'
c      do i=1,N
c        nD = nHtot(i)*Diff(i)
c        if (i==1) then
c          d1 = d1l(1)*xx(1)+d1m(1)*xx(2)+d1r(1)*xx(3)
c        else if (i==N) then
c          d1 = d1l(N)*xx(N-2)+d1m(N)*xx(N-1)+d1r(N)*xx(N)
c        else
c          d1 = d1l(i)*xx(i-1)+d1m(i)*xx(i)+d1r(i)*xx(i+1)
c        endif
c        jdiff = -nD*d1
c        print'(I4,0pF12.5,99(1pE12.4))',i,zz(i)/Hp,xx(i),jdiff
c      enddo
c
c      ntot2 = 0.0
c      do i=1,N-1
c        ntot2 = ntot2 + 0.5*(nHtot(i)*xx(i)+nHtot(i+1)*xx(i+1))
c     >                     *(zz(i+1)-zz(i))
c      enddo
c      print*,ntot,ntot2

      end


************************************************************************
      subroutine DIFFUSION_TINDEPENDENT
************************************************************************
c      use NATURE,ONLY: pi
c      use GRID,ONLY: N=>Npoints,zz,xx,d1l,d1m,d1r,d2l,d2m,d2r
c      use PARAMETERS,ONLY: Hp,bc_low,bc_high,init,tnull,
c     >                     influx,outflux,inrate,outrate,vin,vout
c      use STRUCT,ONLY: Diff,nHtot
c      implicit none
c
c      integer :: i
c      real*8,dimension(N,N) :: A
c      real*8,dimension(N) :: rest
c      real*8 :: nD,d1nD,d1,jdiff
c
c      !-----------------------------------------------------------
c      ! ***  -d/dz(nHtot*Diff*dx/dz) = source terms [1/cm3/s]  ***
c      !-----------------------------------------------------------
c      !-------------------------
c      ! ***    A x = rest    ***
c      !-------------------------
c      A(:,:) = 0.0
c      do i=2,N-1
c        nD   = nHtot(i)*Diff(i)
c        d1nD = d1l(i)*nHtot(i-1)*Diff(i-1)
c     >       + d1m(i)*nHtot(i)  *Diff(i)
c     >       + d1r(i)*nHtot(i+1)*Diff(i+1)
c        A(i,i-1) = A(i,i-1) - d1nD*d1l(i) - nD*d2l(i)
c        A(i,i)   = A(i,i)   - d1nD*d1m(i) - nD*d2m(i)
c        A(i,i+1) = A(i,i+1) - d1nD*d1r(i) - nD*d2r(i)
c      enddo
c
c      !------------------------------
c      ! ***  boundary conditions  ***
c      !------------------------------
c      rest(:) = 0.0
c      nD = nHtot(1)*Diff(1)
c      if (bc_low==1) then
c        A(1,1)  = 1.d0
c        rest(1) = xx(1)
c      else if (bc_low==2) then
c        A(1,1)  = d1l(1)*nD
c        A(1,2)  = d1m(1)*nD
c        A(1,3)  = d1r(1)*nD
c        rest(1) = -influx
c      else if (bc_low==3) then
c        A(1,1)  = d1l(1)*nD + nHtot(1)*inrate*vin
c        A(1,2)  = d1m(1)*nD
c        A(1,3)  = d1r(1)*nD
c        rest(1) = 0.d0
c      endif
c      nD = nHtot(N)*Diff(N)
c      if (bc_high==1) then
c        A(N,N)  = 1.d0
c        rest(N) = xx(N)
c      else if (bc_high==2) then
c        A(N,N-2) = d1l(N)*nD
c        A(N,N-1) = d1m(N)*nD
c        A(N,N)   = d1r(N)*nD
c        rest(N)  = -outflux
c      else if (bc_high==3) then
c        A(N,N-2) = d1l(N)*nD
c        A(N,N-1) = d1m(N)*nD
c        A(N,N)   = d1r(N)*nD + nHtot(N)*outrate*vout
c        rest(N)  = 0.d0
c      endif
c
c      !-------------------------
c      ! *** add source terms ***
c      !-------------------------
c      !rest(N/3) = rest(N/3) + 500.0
c
c      !----------------------------------
c      ! *** solve the matrix equation ***
c      !----------------------------------
c      call GAUSS8(N,N,A,xx,rest)
c
c      open(unit=1,file='out.dat',status='replace')
c      write(1,*) N,init
c      write(1,'(9999(1pE16.8))') zz
c      write(1,'("time[s]=",1pE12.5)') 1.E+40
c      write(1,'(9999(1pE16.8))') (MAX(xx(i),1.E-99),i=1,N)
c      close(1)
c
c      print*
c      print'(A4,99(A12))','grid','zz/Hp','eps','jdiff'
c      do i=1,N
c        nD = nHtot(i)*Diff(i)
c        if (i==1) then
c          d1 = d1l(1)*xx(1)+d1m(1)*xx(2)+d1r(1)*xx(3)
c        else if (i==N) then
c          d1 = d1l(N)*xx(N-2)+d1m(N)*xx(N-1)+d1r(N)*xx(N)
c        else
c          d1 = d1l(i)*xx(i-1)+d1m(i)*xx(i)+d1r(i)*xx(i+1)
c        endif
c        jdiff = -nD*d1
c        print'(I4,0pF12.5,99(1pE12.4))',i,zz(i)/Hp,xx(i),jdiff
c      enddo

      end
