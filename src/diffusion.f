************************************************************************
      subroutine DIFFUSION_EXPLICIT
************************************************************************
      use NATURE,ONLY: pi
      use GRID,ONLY: N=>Npoints,zz,xx,d1l,d1m,d1r,d2l,d2m,d2r
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
        dt = MIN(dt,0.5*(zz(i)-zz(i-1))**2/D)
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
        !-------------------------------------------
        ! ***  d/dt(nH*x) = d/dz(nH*Diff*dx/dz)  ***
        !-------------------------------------------
        rate(:) = 0.0
        do i=2,N-1
          nD   = nHtot(i)*Diff(i)  
          d1   = d1l(i)*xx(i-1) + d1m(i)*xx(i) + d1r(i)*xx(i+1)
          d2   = d2l(i)*xx(i-1) + d2m(i)*xx(i) + d2r(i)*xx(i+1)
          d1nD = d1l(i)*nHtot(i-1)*Diff(i-1) 
     >         + d1m(i)*nHtot(i)  *Diff(i) 
     >         + d1r(i)*nHtot(i+1)*Diff(i+1) 
          rate(i) = nD*d2 + d1nD*d1
        enddo

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
        if (bc_low==1) then
          influx = -nD*(d1l(1)*xx(1) + d1m(1)*xx(2) + d1r(1)*xx(3))
        else if (bc_low==2) then   
          xx(1) = (-influx/nD - d1m(1)*xx(2) - d1r(1)*xx(3))/d1l(1) 
        else if (bc_low==3) then   
          xx(1) = -(d1m(1)*xx(2) + d1r(1)*xx(3))
     >            /(d1l(1) + inrate*vin/Diff(1))
        endif  
        nD = nHtot(N)*Diff(N)  
        if (bc_high==1) then
          outflux = -nD*(d1l(N)*xx(N-2) + d1m(N)*xx(N-1) + d1r(N)*xx(N)) 
        else if (bc_high==2) then   
          xx(N) = (-outflux/nD - d1l(N)*xx(N-2) - d1m(N)*xx(N-1))/d1r(N)
        else if (bc_high==3) then   
          outflux = nHtot(N)*xx(N)*outrate*vout  
          xx(N) = -(d1l(N)*xx(N-2) + d1m(N)*xx(N-1))
     >            /(d1r(N) + outrate*vout/Diff(N))
        endif   
        ntot = ntot + (influx - outflux)*dt
        
        if (time>outtime(Nout)) then
          write(*,'(I8," output t=",1pE11.3," s")') it,time
          write(1,'("time[s]=",1pE12.5)') time 
          write(1,'(9999(1pE16.8))') (MAX(xx(i),1.E-99),i=1,N)
          Nout = Nout + 1
        endif  

        if (time>tend) exit

      enddo  

      print*
      print'(A4,99(A12))','grid','zz/Hp','eps','jdiff'
      do i=1,N
        nD = nHtot(i)*Diff(i)
        if (i==1) then 
          d1 = d1l(1)*xx(1)+d1m(1)*xx(2)+d1r(1)*xx(3) 
        else if (i==N) then 
          d1 = d1l(N)*xx(N-2)+d1m(N)*xx(N-1)+d1r(N)*xx(N)
        else
          d1 = d1l(i)*xx(i-1)+d1m(i)*xx(i)+d1r(i)*xx(i+1)
        endif   
        jdiff = -nD*d1
        print'(I4,0pF12.5,99(1pE12.4))',i,zz(i)/Hp,xx(i),jdiff
      enddo

      ntot2 = 0.0
      do i=1,N-1
        ntot2 = ntot2 + 0.5*(nHtot(i)*xx(i)+nHtot(i+1)*xx(i+1))
     >                     *(zz(i+1)-zz(i))
      enddo  
      print*,ntot,ntot2

      end


************************************************************************
      subroutine DIFFUSION_IMPLICIT
************************************************************************
      use NATURE,ONLY: pi
      use GRID,ONLY: N=>Npoints,zz,xx,d1l,d1m,d1r,d2l,d2m,d2r
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
        dt = MIN(dt,0.5*(zz(i)-zz(i-1))**2/D)
      enddo
      print*,"explicit timestep=",dt
      dt = dt*tfac
      dt = dt*(1.d0+1.d-12)
      print*," applied timestep=",dt

      !-----------------------------
      ! ***  fill in big matrix  ***
      !-----------------------------
      A(:,:) = 0.d0
      do i=2,N-1
        nD   = nHtot(i)*Diff(i)  
        d1nD = d1l(i)*nHtot(i-1)*Diff(i-1) 
     >       + d1m(i)*nHtot(i)  *Diff(i) 
     >       + d1r(i)*nHtot(i+1)*Diff(i+1) 
        A(i,i-1) = A(i,i-1) - dt*( d1nD*d1l(i) + nD*d2l(i) )
        A(i,i)   = A(i,i)   - dt*( d1nD*d1m(i) + nD*d2m(i) )
        A(i,i+1) = A(i,i+1) - dt*( d1nD*d1r(i) + nD*d2r(i) )
      enddo
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
      if (bc_low==1) then
        !--- nothing to do 
      else if (bc_low==2) then
        A(1,1) = 1.d0
        A(1,2) = d1m(1)/d1l(1)
        A(1,3) = d1r(1)/d1l(1)
      else if (bc_low==3) then
        A(1,1) = 1.d0+inrate*vin/Diff(1)/d1l(1) 
        A(1,2) = d1m(1)/d1l(1)
        A(1,3) = d1r(1)/d1l(1)
      endif  
      if (bc_high==1) then
        !--- nothing to do 
      else if (bc_high==2) then   
        A(N,N-2) = d1l(N)/d1r(N)
        A(N,N-1) = d1m(N)/d1r(N)
        A(N,N)   = 1.d0
      else if (bc_high==3) then   
        A(N,N-2) = d1l(N)/d1r(N)
        A(N,N-1) = d1m(N)/d1r(N)
        A(N,N)   = 1.d0+outrate*vout/Diff(N)/d1r(N) 
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

      ntot = 0.0
      do i=1,N-1
        ntot = ntot + 0.5*(nHtot(i)*xx(i)+nHtot(i+1)*xx(i+1))
     >                   *(zz(i+1)-zz(i))
      enddo  

      time = tnull
      tend = outtime(Nout)
      Nout  = 1
      print*
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
        rest = xx
        nD = nHtot(1)*Diff(1)  
        if (bc_low==1) then
          influx  = -nD*(d1l(1)*xx(1) + d1m(1)*xx(2) + d1r(1)*xx(3))
        else if (bc_low==2) then
          rest(1) = -influx/nD/d1l(1) 
        else if (bc_low==3) then
          influx  = inrate*nHtot(1)*xx(1)*vin 
          rest(1) = 0.d0
        endif
        nD = nHtot(N)*Diff(N)  
        if (bc_high==1) then
          outflux = -nD*(d1l(N)*xx(N-2) + d1m(N)*xx(N-1) + d1r(N)*xx(N)) 
        else if (bc_high==2) then
          rest(N) = -outflux/nD/d1r(N)           
        else if (bc_high==3) then
          outflux = outrate*nHtot(N)*xx(N)*vout 
          rest(N) = 0.d0
        endif
        !ntot = ntot + (influx-outflux)*dt  

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
          influx  = -nD*(d1l(1)*xx(1) + d1m(1)*xx(2) + d1r(1)*xx(3))
        else if (bc_low==3) then  
          influx  = inrate*nHtot(1)*xx(1)*vin 
        endif
        nD = nHtot(N)*Diff(N)  
        if (bc_high==1) then
          outflux = -nD*(d1l(N)*xx(N-2) + d1m(N)*xx(N-1) + d1r(N)*xx(N)) 
        else if (bc_high==3) then  
          outflux = outrate*nHtot(N)*xx(N)*vout 
        endif
        ntot = ntot + (influx-outflux)*dt  
        time = time + dt
        write(*,'(TL10,I8,A,$)') it,CR

        if (time>outtime(Nout)) then
          write(*,'(I8," output t=",1pE11.3," s")') it,time
          write(1,'("time[s]=",1pE12.5)') time 
          write(1,'(9999(1pE16.8))') (MAX(xx(i),1.E-99),i=1,N)
          Nout = Nout + 1
        endif  

        if (time>tend) exit

      enddo  

      print*
      print'(A4,99(A12))','grid','zz/Hp','eps','jdiff'
      do i=1,N
        nD = nHtot(i)*Diff(i)
        if (i==1) then 
          d1 = d1l(1)*xx(1)+d1m(1)*xx(2)+d1r(1)*xx(3) 
        else if (i==N) then 
          d1 = d1l(N)*xx(N-2)+d1m(N)*xx(N-1)+d1r(N)*xx(N)
        else
          d1 = d1l(i)*xx(i-1)+d1m(i)*xx(i)+d1r(i)*xx(i+1)
        endif   
        jdiff = -nD*d1
        print'(I4,0pF12.5,99(1pE12.4))',i,zz(i)/Hp,xx(i),jdiff
      enddo

      ntot2 = 0.0
      do i=1,N-1
        ntot2 = ntot2 + 0.5*(nHtot(i)*xx(i)+nHtot(i+1)*xx(i+1))
     >                     *(zz(i+1)-zz(i))
      enddo  
      print*,ntot,ntot2

      end


************************************************************************
      subroutine DIFFUSION_TINDEPENDENT
************************************************************************
      use NATURE,ONLY: pi
      use GRID,ONLY: N=>Npoints,zz,xx,d1l,d1m,d1r,d2l,d2m,d2r
      use PARAMETERS,ONLY: Hp,bc_low,bc_high,init,tnull,
     >                     influx,outflux,inrate,outrate,vin,vout
      use STRUCT,ONLY: Diff,nHtot
      implicit none

      integer :: i
      real*8,dimension(N,N) :: A
      real*8,dimension(N) :: rest
      real*8 :: nD,d1nD,d1,jdiff

      !-----------------------------------------------------------
      ! ***  -d/dz(nHtot*Diff*dx/dz) = source terms [1/cm3/s]  ***
      !-----------------------------------------------------------
      !-------------------------
      ! ***    A x = rest    ***
      !-------------------------
      A(:,:) = 0.0
      do i=2,N-1
        nD   = nHtot(i)*Diff(i)  
        d1nD = d1l(i)*nHtot(i-1)*Diff(i-1) 
     >       + d1m(i)*nHtot(i)  *Diff(i) 
     >       + d1r(i)*nHtot(i+1)*Diff(i+1) 
        A(i,i-1) = A(i,i-1) - d1nD*d1l(i) - nD*d2l(i)
        A(i,i)   = A(i,i)   - d1nD*d1m(i) - nD*d2m(i)
        A(i,i+1) = A(i,i+1) - d1nD*d1r(i) - nD*d2r(i)
      enddo

      !---------------------------------
      ! *** (2) fill in source terms ***
      !---------------------------------
      rest(:) = 0.0

      !------------------------------
      ! ***  boundary conditions  ***
      !------------------------------
      nD = nHtot(1)*Diff(1)  
      if (bc_low==1) then
        A(1,1)  = 1.d0
        rest(1) = xx(1) 
      else if (bc_low==2) then
        A(1,1)  = d1l(1)*nD
        A(1,2)  = d1m(1)*nD
        A(1,3)  = d1r(1)*nD
        rest(1) = -influx
      else if (bc_low==3) then
        A(1,1)  = d1l(1)*nD + nHtot(1)*inrate*vin
        A(1,2)  = d1m(1)*nD
        A(1,3)  = d1r(1)*nD
        rest(1) = 0.d0
      endif  
      nD = nHtot(N)*Diff(N)  
      if (bc_high==1) then
        A(N,N)  = 1.d0
        rest(N) = xx(N) 
      else if (bc_high==2) then   
        A(N,N-2) = d1l(N)*nD
        A(N,N-1) = d1m(N)*nD
        A(N,N)   = d1r(N)*nD
        rest(N)  = -outflux
      else if (bc_high==3) then   
        A(N,N-2) = d1l(N)*nD
        A(N,N-1) = d1m(N)*nD
        A(N,N)   = d1r(N)*nD + nHtot(N)*outrate*vout
        rest(N)  = 0.d0
      endif   

      !----------------------------------
      ! *** solve the matrix equation ***
      !----------------------------------
      call GAUSS8(N,N,A,xx,rest)

      open(unit=1,file='out.dat',status='replace')
      write(1,*) N,init 
      write(1,'(9999(1pE16.8))') zz
      write(1,'("time[s]=",1pE12.5)') 1.E+40
      write(1,'(9999(1pE16.8))') (MAX(xx(i),1.E-99),i=1,N)
      close(1)

      print*
      print'(A4,99(A12))','grid','zz/Hp','eps','jdiff'
      do i=1,N
        nD = nHtot(i)*Diff(i)
        if (i==1) then 
          d1 = d1l(1)*xx(1)+d1m(1)*xx(2)+d1r(1)*xx(3) 
        else if (i==N) then 
          d1 = d1l(N)*xx(N-2)+d1m(N)*xx(N-1)+d1r(N)*xx(N)
        else
          d1 = d1l(i)*xx(i-1)+d1m(i)*xx(i)+d1r(i)*xx(i+1)
        endif   
        jdiff = -nD*d1
        print'(I4,0pF12.5,99(1pE12.4))',i,zz(i)/Hp,xx(i),jdiff
      enddo

      end
