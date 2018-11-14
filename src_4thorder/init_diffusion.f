************************************************************************
      subroutine INIT_BMATRIX(N,bc_low,bc_high,dt,B)
************************************************************************
      use GRID,ONLY: zz,d1l1,d1l2,d1m,d1r1,d1r2,
     >                  d2l1,d2l2,d2m,d2r1,d2r2
      use STRUCT,ONLY: Diff,nHtot
      use PARAMETERS,ONLY: influx,outflux,inrate,outrate,vin,vout
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(IN) :: N,bc_low,bc_high
      real*8,intent(IN)  :: dt
      real*8,intent(OUT) :: B(N,N)
      real*8 :: nD,d1nD,err,A(N,N),sum(N,N)
      real(kind=qp),dimension(N,N) :: Awork
      real(kind=qp) :: det(2),work(N)
      integer :: i,j,k,ipvt(N),info
      logical :: check=.false.

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
      nD   = nHtot(2)*Diff(2)
      d1nD = d1l2(2)*nHtot(1)*Diff(1)
     >     + d1l1(2)*nHtot(2)*Diff(2)
     >     + d1m(2) *nHtot(3)*Diff(3)
     >     + d1r1(2)*nHtot(4)*Diff(4)
     >     + d1r2(2)*nHtot(5)*Diff(5)
      A(2,1) = A(2,1) - dt*( d1nD*d1l2(2) + nD*d2l2(2) )
      A(2,2) = A(2,2) - dt*( d1nD*d1l1(2) + nD*d2l1(2) )
      A(2,3) = A(2,3) - dt*( d1nD*d1m(2)  + nD*d2m(2)  )
      A(2,4) = A(2,4) - dt*( d1nD*d1r1(2) + nD*d2r1(2) )
      A(2,5) = A(2,5) - dt*( d1nD*d1r2(2) + nD*d2r2(2) )
      !--- i=N-1 ---
      nD   = nHtot(N-1)*Diff(N-1)
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
      !--- i=1,N come later as boundary conditions ---
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
      if (bc_low==1) then                   !fixed concentration
        !--- nothing to do
      else if (bc_low==2) then              !fixed influx
        A(1,1) = 1.d0
        A(1,2) = d1l1(1)/d1l2(1)
        A(1,3) =  d1m(1)/d1l2(1)
        A(1,4) = d1r1(1)/d1l2(1)
        A(1,5) = d1r2(1)/d1l2(1)
      else if (bc_low==3) then              !fixed inflow rate
        A(1,1) = 1.d0+inrate*vin/Diff(1)/d1l2(1)
        A(1,2) = d1l1(1)/d1l2(1)
        A(1,3) =  d1m(1)/d1l2(1)
        A(1,4) = d1r1(1)/d1l2(1)
        A(1,5) = d1r2(1)/d1l2(1)
      endif
      if (bc_high==1) then                  !fixed concentration
        !--- nothing to do
      else if (bc_high==2) then             !fixed outflux
        A(N,N-4) = d1l2(N)/d1r2(N)
        A(N,N-3) = d1l1(N)/d1r2(N)
        A(N,N-2) =  d1m(N)/d1r2(N)
        A(N,N-1) = d1r1(N)/d1r2(N)
        A(N,N)   = 1.d0
      else if (bc_high==3) then             !fixed outflow rare
        A(N,N-4) = d1l2(N)/d1r2(N)
        A(N,N-3) = d1l1(N)/d1r2(N)
        A(N,N-2) =  d1m(N)/d1r2(N)
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
      
      end
