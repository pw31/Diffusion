************************************************************************
      subroutine DIFFUSION_EXPLICIT
************************************************************************
      use NATURE,ONLY: pi
      use GRID,ONLY: N=>Npoints,zz,zweight,xx,d1l2,d1l1,d1m,d1r1,
     >d1r2,d2l2,d2l1,d2m,d2r1,d2r2
      use PARAMETERS,ONLY: Hp,tnull,init,bc_low,bc_high,
     >                     influx,outflux,inrate,outrate,vin,vout,
     >                     Nout,outtime
      use STRUCT,ONLY: Diff,nHtot
      implicit none
      integer :: i,it
      real*8,dimension(N) :: rate,xold
      real*8 :: D,nD,d1,d2,d1nD,jdiff,time,dt,tend
      real*8 :: ww,AA,z0,x1ana,xNana
      real*8 :: ntot,ntot2,dNcol,influx0
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
      print*,"explicit diffusion timestep[s]=",dt
      print*
      
      ntot = 0.d0
      do i=1,N
        ntot = ntot + zweight(i)*nHtot(i)*xx(i)
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
        else if (bc_low==2) then                 ! fixed influx
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
        else if (bc_high==2) then                ! fixed outflux
          xx(N) = (-outflux/nD -d1l2(N)*xx(N-4)
     >                         -d1l1(N)*xx(N-3) 
     >                         -d1m(N) *xx(N-2)
     >                         -d1r1(N)*xx(N-1) )/d1r2(N)
        else if (bc_high==3) then                ! fixed outrate
          outflux = nHtot(N)*xx(N)*outrate*vout  
          xx(N) = -( d1l2(N)*xx(N-4)
     >              +d1l1(N)*xx(N-3)
     >              +d1m(N) *xx(N-2) 
     >              +d1r1(N)*xx(N-1) )/(d1r2(N) + outrate*vout/Diff(N))
        endif   
        ntot = ntot + (influx - outflux)*0.5*dt
        influx0 = influx
        
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
        xold = xx
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

        !dNcol = 0.d0
        !do i=1,N
        !  dNcol = dNcol + zweight(i)*nHtot(i)*(xx(i)-xold(i))
        !enddo  
        !print*,dt,bc_low,bc_high
        !print'(7(1pE12.5))',xold(1:7)
        !print'(7(1pE12.5))',xx(1:7)
        !print*,0.5*(influx+influx0),dNcol/dt

        if (.true..or.time>outtime(Nout)) then
          write(*,'(I8," output t=",1pE14.6," s")') it,time
          write(1,'("time[s]=",1pE12.5)') time 
          write(1,'(9999(1pE16.8))') (MAX(xx(i),1.E-99),i=1,N)
          Nout = Nout + 1
          ntot2 = 0.0
          do i=1,N
            ntot2 = ntot2 + zweight(i)*nHtot(i)*xx(i)
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
      use NATURE,ONLY: pi
      use GRID,ONLY: N=>Npoints,zz,xx,d1l2,d1l1,d1m,d1r1,d1r2,
     >              d2l2,d2l1,d2m,d2r1,d2r2,zweight
      use PARAMETERS,ONLY: Hp,bc_low,bc_high,init,tnull,
     >                     influx,outflux,inrate,outrate,vin,vout,
     >                     Nout,outtime,tfac
      use STRUCT,ONLY: Diff,nHtot
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer :: i,j,k,it,ipvt(N),info
      real*8,dimension(N,N) :: A,B,sum
      real(kind=qp),dimension(N,N) :: Awork
      real(kind=qp) :: det(2),work(N)
      real*8,dimension(N) :: xnew,rest
      real*8 :: D,nD,d1,d1nD,jdiff,time,dt,tend
      real*8 :: z0,ww,AA,x1ana,xNana,ntot,ntot2,err
      logical :: check=.true.
      character :: CR = CHAR(13)

      !-------------------
      ! ***  timestep  ***
      !-------------------
      dt = 9.D+99
      do i=2,N
        D  = 0.5*(Diff(i-1)+Diff(i))
        dt = MIN(dt,0.33*(zz(i)-zz(i-1))**2/D)
      enddo
      print*,"explicit timestep=",dt
      dt = dt*tfac
      dt = dt*(1.d0+1.d-12)
      print*," applied timestep=",dt
c
      !-----------------------------
      ! ***  fill in big matrix  ***
      !-----------------------------
      A(:,:) = 0.d0
      do i=3,N-2
        nD   = nHtot(i)*Diff(i)
        d1nD = d1l2(i)*nHtot(i-2)*Diff(i-2)
     >       + d1l1(i)*nHtot(i-1)*Diff(i-1)
     >       + d1m(i) *nHtot(i)  *Diff(i)
     >       + d1r1(i)*nHtot(i+1)*Diff(i+1)
     >       + d1r2(i)*nHtot(i+2)*Diff(i+2)
        A(i,i-2) = A(i,i-2) - dt*( d1nD*d1l2(i) + nD*d2l2(i) )
        A(i,i-1) = A(i,i-1) - dt*( d1nD*d1l1(i) + nD*d2l1(i) )
        A(i,i)   = A(i,i)   - dt*( d1nD*d1m(i)  + nD*d2m(i) )
        A(i,i+1) = A(i,i+1) - dt*( d1nD*d1r1(i) + nD*d2r1(i) )
        A(i,i+2) = A(i,i+2) - dt*( d1nD*d1r2(i) + nD*d2r2(i) )
      enddo

      !--- i=2 ---
      d1nD = d1l2(2)*nHtot(1)*Diff(1)
     >       + d1l1(2)*nHtot(2)*Diff(2)
     >       + d1m(2) *nHtot(3)*Diff(3)
     >       + d1r1(2)*nHtot(4)*Diff(4)
     >       + d1r2(2)*nHtot(5)*Diff(5)
      A(2,1) = A(2,1) - dt*( d1nD*d1l2(2) + nD*d2l2(2) )
      A(2,2) = A(2,2) - dt*( d1nD*d1l1(2) + nD*d2l1(2) )
      A(2,3) = A(2,3) - dt*( d1nD*d1m(2)  + nD*d2m(2)  )
      A(2,4) = A(2,4) - dt*( d1nD*d1r1(2) + nD*d2r1(2) )
      A(2,5) = A(2,5) - dt*( d1nD*d1r2(2) + nD*d2r2(2) )
      !--- i=N-1 ---
      d1nD = d1l2(N-1)*nHtot(N-4)*Diff(N-4)
     >     + d1l1(N-1)*nHtot(N-3)*Diff(N-3)
     >     + d1m(N-1) *nHtot(N-2)*Diff(N-2)
     >     + d1r1(N-1)*nHtot(N-1)*Diff(N-1)
     >     + d1r2(N-1)*nHtot(N)  *Diff(N)
      A(N-1,N-4) = A(N-1,N-4) - dt*( d1nD*d1l2(N-1) + nD*d2l2(N-1) )
      A(N-1,N-3) = A(N-1,N-3) - dt*( d1nD*d1l1(N-1) + nD*d2l1(N-1) )
      A(N-1,N-2) = A(N-1,N-2) - dt*( d1nD*d1m(N-1)  + nD*d2m(N-1)  )
      A(N-1,N-1) = A(N-1,N-1) - dt*( d1nD*d1r1(N-1) + nD*d2r1(N-1) )
      A(N-1,N)   = A(N-1,N)   - dt*( d1nD*d1r2(N-1) + nD*d2r2(N-1) )

c no entry yet for i = 1,N
      do i=1,N
        A(i,:) = A(i,:)/nHtot(i)    ! unitless
      enddo
      !--------------------------
      ! ***  add unit matrix  ***
      !--------------------------
      do i=1,N
        A(i,i) = A(i,i) + 1.d0
      enddo
      !------------------------------
      ! ***  boundary conditions  ***
      !------------------------------
      if (bc_low==1) then                   !fiexed concentration
        !--- nothing to do
      else if (bc_low==2) then              !fixed flux
        A(1,1) = 1.d0
        A(1,2) = d1l1(1)/d1l2(1)
        A(1,3) = d1m(1)/d1l2(1)
        A(1,4) = d1r1(1)/d1l2(1)
        A(1,5) = d1r2(1)/d1l2(1)
      else if (bc_low==3) then              !fixed outflow rate
        A(1,1) = 1.d0+inrate*vin/Diff(1)/d1l2(1)
        A(1,2) = d1l1(1)/d1l2(1)
        A(1,3) = d1m(1)/d1l2(1)
        A(1,4) = d1r1(1)/d1l2(1)
        A(1,5) = d1r2(1)/d1l2(1)
      endif
      if (bc_high==1) then                  !fixed concentration
        !--- nothing to do
      else if (bc_high==2) then             !fixed flux
        A(N,N-4) = d1l2(N)/d1r2(N)
        A(N,N-4) = d1l1(N)/d1r2(N)
        A(N,N-2) = d1m(N)/d1r2(N)
        A(N,N-1) = d1r1(N)/d1r2(N)
        A(N,N)   = 1.d0
      else if (bc_high==3) then             !fixed outflow rare
        A(N,N-4) = d1l2(N)/d1r2(N)
        A(N,N-4) = d1l1(N)/d1r2(N)
        A(N,N-2) = d1m(N)/d1r2(N)
        A(N,N-1) = d1r1(N)/d1r2(N)
        A(N,N)   = 1.d0+outrate*vout/Diff(N)/d1r2(N)
      endif

      !------------------------
      ! ***  invert matrix  ***
      !------------------------
      Awork = A
      call QGEFA ( Awork, N, N, ipvt, info )
      call QGEDI ( Awork, N, N, ipvt, det, work, 1 )
      B = Awork
      if (info.ne.0) then
        print*,"*** singular matrix in QGEFA: info=",info
        stop
      endif
      if (check) then
        do i=1,N
          write(99,'(9999(1pE11.3))') (A(i,j),j=1,N)
        enddo
        write(99,*)
        do i=1,N
          write(99,'(9999(1pE11.3))') (B(i,j),j=1,N)
        enddo
        !--- test A*B=1 ---
        err = 0.d0
        write(99,*)
        do i=1,N
          do j=1,N
            sum(i,j) = 0.d0
            do k=1,N
              sum(i,j) = sum(i,j) + A(i,k)*B(k,j)
            enddo
            if (i==j) then
              err = max(err,ABS(sum(i,j)-1.d0))
            else
              err = max(err,ABS(sum(i,j)))
            endif
          enddo
          write(99,'(9999(1pE11.3))') (sum(i,j),j=1,N)
        enddo
        write(99,*) "maximum error=",err
        print*,"matrix inversion error=",err
      endif
c
      ntot = 0.0
      do i=1,N-1
        ntot = ntot + 0.5*(nHtot(i)*xx(i)+nHtot(i+1)*xx(i+1))
     >                   *(zz(i+1)-zz(i))
      enddo
c
      time = tnull
      tend = outtime(Nout)
      Nout  = 1
      print*
      open(unit=1,file='out.dat',status='replace')
      write(1,*) N,init
      write(1,'(9999(1pE16.8))') zz
      write(1,'("time[s]=",1pE12.5)') time
      write(1,'(9999(1pE16.8))') (MAX(xx(i),1.E-99),i=1,N)
c
      do it=1,9999999
c
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
        rest = xx
        nD = nHtot(1)*Diff(1)
        if (bc_low==1) then                     !fixed concentration
          influx  = -nD*(d1l2(1)*xx(1)
     >                 + d1l1(1)*xx(2)
     >                 + d1m (1)*xx(3)
     >                 + d1r1(1)*xx(4)
     >                 + d1r2(1)*xx(5))
        else if (bc_low==2) then                !fixed flux
          rest(1) = -influx/nD/d1l2(1)
        else if (bc_low==3) then                !fixed outflow rate
          influx  = inrate*nHtot(1)*xx(1)*vin
          rest(1) = 0.d0
        endif
        nD = nHtot(N)*Diff(N)
        if (bc_high==1) then                    !fixed concentration
          outflux = -nD*(d1l2(N)*xx(N-4)
     >                 + d1l1(N)*xx(N-3)
     >                 + d1m(N) *xx(N-2)
     >                 + d1r1(N)*xx(N-1)
     >                 + d1r2(N)*xx(N))
        else if (bc_high==2) then               !fixed flux
          rest(N) = -outflux/nD/d1r2(N)
        else if (bc_high==3) then               !fixed outflow rate
          outflux = outrate*nHtot(N)*xx(N)*vout
          rest(N) = 0.d0
        endif
        !ntot = ntot + (influx-outflux)*dt
c
        !-----------------------
        ! ***  source terms  ***
        !-----------------------
        !rest(N/3) = rest(N/3) + 500.0*dt
c
        !----------------------------
        ! ***  implicit timestep  ***
        !----------------------------
        do i=1,N
          xnew(i) = 0.d0
          do j=1,N
            xnew(i) = xnew(i) + B(i,j)*rest(j)
          enddo
        enddo
        xx(:) = xnew(:)
        nD = nHtot(1)*Diff(1)
        if (bc_low==1) then
          influx  = -nD*(d1l2(1)*xx(1)
     >                 + d1l1(1)*xx(2)
     >                 + d1m(1) *xx(3)
     >                 + d1r1(1)*xx(4)
     >                 + d1r2(1)*xx(5))
        else if (bc_low==3) then
          influx  = inrate*nHtot(1)*xx(1)*vin
        endif
        nD = nHtot(N)*Diff(N)
        if (bc_high==1) then
          outflux = -nD*(d1l2(N)*xx(N-4)
     >                 + d1l1(N)*xx(N-3)
     >                 + d1m(N) *xx(N-2)
     >                 + d1r1(N)*xx(N-1)
     >                 + d1r2(N)*xx(N))
        else if (bc_high==3) then
          outflux = outrate*nHtot(N)*xx(N)*vout
        endif
        ntot = ntot + (influx-outflux)*dt
        time = time + dt
        write(*,'(TL10,I8,A,$)') it,CR
c
        !if (time>outtime(Nout)) then
        !  write(*,'(I8," output t=",1pE11.3," s")') it,time
        !  write(1,'("time[s]=",1pE12.5)') time
        !  write(1,'(9999(1pE16.8))') (MAX(xx(i),1.E-99),i=1,N)
        !  Nout = Nout + 1
        !endif
        if (time>outtime(Nout)) then
            write(*,'(I8," output t=",1pE14.6," s")') it,time
            write(1,'("time[s]=",1pE12.5)') time
            write(1,'(9999(1pE16.8))') (MAX(xx(i),1.E-99),i=1,N)
            Nout = Nout + 1
            ntot2 = 0.0
            do i=1,N
                ntot2 = ntot2 + zweight(i)*nHtot(i)*xx(i)
            enddo
            print'("  total=",2(1pE14.6)," , dev=",0pF8.5,"%")',
     >         ntot,ntot2,(ntot/ntot2-1.0)*100.0
        endif
c
        if (time>tend) exit

      enddo

      print*
      print'(A4,99(A12))','grid','zz/Hp','eps','jdiff'
      do i=1,N
        nD = nHtot(i)*Diff(i)
        if (i==1) then
          d1 = d1l2(1)*xx(1)
     >       + d1l1(1)*xx(2)
     >       + d1m(1) *xx(3)
     >       + d1r1(1)*xx(4)
     >       + d1r2(1)*xx(5)
        else if (i==2) then
          d1 = d1l2(2)*xx(1)
     >       + d1l1(2)*xx(2)
     >       + d1m(2)*xx(3)
     >       + d1r1(2)*xx(4)
     >       + d1r2(2)*xx(5)
        else if (i==N-1) then
          d1 = d1l1(N-1)*xx(N-4)
     >       + d1l1(N-1)*xx(N-3)
     >       + d1m(N-1) *xx(N-2)
     >       + d1r1(N-1)*xx(N-1)
     >       + d1r2(N-1)*xx(N)
        else if (i==N) then
          d1 = d1l1(N)*xx(N-4)
     >       + d1l1(N)*xx(N-3)
     >       + d1m(N) *xx(N-2)
     >       + d1r1(N)*xx(N-1)
     >       + d1r2(N)*xx(N)
        else
          d1 = d1l2(i)*xx(i-2)
     >       + d1l1(i)*xx(i-1)
     >       + d1m(i) *xx(i)
     >       + d1r1(i)*xx(i+1)
     >       + d1r2(i)*xx(i+2)
        endif
        jdiff = -nD*d1
        print'(I4,0pF12.5,99(1pE12.4))',i,zz(i)/Hp,xx(i),jdiff
      enddo
c
      !ntot2 = 0.0
      !do i=1,N-1
      !  ntot2 = ntot2 + 0.5*(nHtot(i)*xx(i)+nHtot(i+1)*xx(i+1))
      !>                     *(zz(i+1)-zz(i))
      !enddo
      !print*,ntot,ntot2
      print'("  total=",2(1pE14.6)," , dev=",0pF8.5,"%")',
     >         ntot,ntot2,(ntot/ntot2-1.0)*100.0
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
