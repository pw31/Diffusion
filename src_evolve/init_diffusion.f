************************************************************************
      subroutine INIT_DIFFUSION(N,bc_low,dt,BB)
************************************************************************
      use GRID,ONLY: zz,d1l,d1m,d1r,d2l,d2m,d2r
      use STRUCT,ONLY: nHtot,Diff
      use PARAMETERS,ONLY: bc_high,inrate,outrate,vin,vout
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(in) :: N,bc_low
      real*8,intent(in) :: dt
      real*8,intent(out) :: BB(N,N)
      integer :: i,j,k,ipvt(N),info
      real(kind=qp),dimension(N,N) :: Awork
      real(kind=qp) :: det(2),work(N)
      real*8 :: AA(N,N),sum(N,N)
      real*8 :: df,h1,h2
      real*8 :: D,nD,d1nD,err
      real*8,dimension(N) :: f0,f1,f2
      logical :: check=.false.

      !-----------------------------
      ! ***  fill in big matrix  ***
      !-----------------------------
      AA(:,:) = 0.d0
      do i=2,N-1
        nD   = nHtot(i)*Diff(i)  
        d1nD = d1l(i)*nHtot(i-1)*Diff(i-1) 
     >       + d1m(i)*nHtot(i)  *Diff(i) 
     >       + d1r(i)*nHtot(i+1)*Diff(i+1) 
        AA(i,i-1) = AA(i,i-1) - dt*( d1nD*d1l(i) + nD*d2l(i) )
        AA(i,i)   = AA(i,i)   - dt*( d1nD*d1m(i) + nD*d2m(i) )
        AA(i,i+1) = AA(i,i+1) - dt*( d1nD*d1r(i) + nD*d2r(i) )
      enddo
      do i=1,N
        AA(i,:) = AA(i,:)/nHtot(i)    ! unitless
      enddo  
      !--------------------------
      ! ***  add unit matrix  ***
      !--------------------------
      do i=1,N
        AA(i,i) = AA(i,i) + 1.d0
      enddo 
      !------------------------------
      ! ***  boundary conditions  ***
      !------------------------------
      if (bc_low==1) then
        !--- nothing to do 
      else if (bc_low==2) then
        AA(1,1) = 1.d0
        AA(1,2) = d1m(1)/d1l(1)
        AA(1,3) = d1r(1)/d1l(1)
      else if (bc_low==3) then
        AA(1,1) = 1.d0+inrate*vin/Diff(1)/d1l(1) 
        AA(1,2) = d1m(1)/d1l(1)
        AA(1,3) = d1r(1)/d1l(1)
      endif  
      if (bc_high==1) then
        !--- nothing to do 
      else if (bc_high==2) then   
        AA(N,N-2) = d1l(N)/d1r(N)
        AA(N,N-1) = d1m(N)/d1r(N)
        AA(N,N)   = 1.d0
      else if (bc_high==3) then   
        AA(N,N-2) = d1l(N)/d1r(N)
        AA(N,N-1) = d1m(N)/d1r(N)
        AA(N,N)   = 1.d0+outrate*vout/Diff(N)/d1r(N) 
      endif   

      !------------------------
      ! ***  invert matrix  ***
      !------------------------
      Awork = AA
      call QGEFA ( Awork, N, N, ipvt, info )
      call QGEDI ( Awork, N, N, ipvt, det, work, 1 )
      BB = Awork
      if (info.ne.0) then
        print*,"*** singular matrix in QGEFA: info=",info
        stop
      endif 
  
      if (check) then
        do i=1,N
          write(99,'(9999(1pE11.3))') (AA(i,j),j=1,N)
        enddo
        write(99,*)
        do i=1,N
          write(99,'(9999(1pE11.3))') (BB(i,j),j=1,N)
        enddo
        !--- test A*B=1 ---
        err = 0.d0
        write(99,*)
        do i=1,N
          do j=1,N
            sum(i,j) = 0.d0
            do k=1,N
              sum(i,j) = sum(i,j) + AA(i,k)*BB(k,j)
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
