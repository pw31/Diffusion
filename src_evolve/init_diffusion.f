************************************************************************
      subroutine INIT_DIFFUSION(el,N,bc_low,bc_high,dt,BB)
************************************************************************
      use GRID,ONLY: zz,d1l,d1m,d1r,d2l,d2m,d2r,dd1l,dd1m,dd1r
      use STRUCT,ONLY: nHtot,Diff
      use JEANS_ESCAPE, ONLY: jpern
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(in) :: el,N,bc_low,bc_high
      real*8,intent(in) :: dt
      real*8,intent(out) :: BB(N+1,N+1)
      integer :: i,j,k,ipvt(N+1),info
      real(kind=qp),dimension(N+1,N+1) :: Awork
      real(kind=qp) :: det(2),work(N+1)
      real*8 :: AA(N+1,N+1),sum(N+1,N+1)
      real*8 :: df,h1,h2
      real*8 :: D,nD,d1nD,err
      real*8,dimension(N+1) :: f0,f1,f2
      logical :: check=.false.

      !-----------------------------
      ! ***  fill in big matrix  ***
      !-----------------------------
      AA(:,:) = 0.d0
      do i=1,N-1
        nD   = nHtot(i)*Diff(i)  
        d1nD = d1l(i)*nHtot(i-1)*Diff(i-1) 
     >       + d1m(i)*nHtot(i)  *Diff(i) 
     >       + d1r(i)*nHtot(i+1)*Diff(i+1) 
        AA(i+1,i+0) = AA(i+1,i+0) - dt*( d1nD*d1l(i) + nD*d2l(i) )
        AA(i+1,i+1) = AA(i+1,i+1) - dt*( d1nD*d1m(i) + nD*d2m(i) )
        AA(i+1,i+2) = AA(i+1,i+2) - dt*( d1nD*d1r(i) + nD*d2r(i) )
      enddo
      do i=0,N-1
        AA(i+1,:) = AA(i+1,:)/nHtot(i)    ! unitless
      enddo  
      !--------------------------
      ! ***  add unit matrix  ***
      !--------------------------
      do i=1,N+1
        AA(i,i) = AA(i,i) + 1.d0
      enddo 
      !------------------------------
      ! ***  boundary conditions  ***
      !------------------------------
      if (bc_low==1) then
        !--- nothing to do 
      else if (bc_low==2) then
        AA(1,1) = 1.d0
        AA(1,2) = dd1m(0)/dd1l(0)
        AA(1,3) = dd1r(0)/dd1l(0)
      endif  
      if (bc_high==1) then
        !--- nothing to do 
      else if (bc_high==2) then   
        AA(N+1,N-1) = d1l(N)/d1r(N)
        AA(N+1,N+0) = d1m(N)/d1r(N)
        AA(N+1,N+1) = 1.d0
      else if (bc_high==3) then   
        AA(N+1,N-1) = d1l(N)/d1r(N)
        AA(N+1,N+0) = d1m(N)/d1r(N)
        AA(N+1,N+1) = 1.d0+jpern(el)/Diff(N)/d1r(N) 
      endif   

      !------------------------
      ! ***  invert matrix  ***
      !------------------------
      Awork = AA
      call QGEFA ( Awork, N+1, N+1, ipvt, info )
      call QGEDI ( Awork, N+1, N+1, ipvt, det, work, 1 )
      BB = Awork
      if (info.ne.0) then
        print*,"*** singular matrix in QGEFA: info=",info
        stop
      endif 
  
      if (check) then
        do i=1,N+1
          write(99,'(9999(1pE11.3))') (AA(i,j),j=1,N+1)
        enddo
        write(99,*)
        do i=1,N+1
          write(99,'(9999(1pE11.3))') (BB(i,j),j=1,N+1)
        enddo
        !--- test A*B=1 ---
        err = 0.d0
        write(99,*)
        do i=1,N+1
          do j=1,N+1
            sum(i,j) = 0.d0
            do k=1,N+1
              sum(i,j) = sum(i,j) + AA(i,k)*BB(k,j)
            enddo 
            if (i==j) then
              err = max(err,ABS(sum(i,j)-1.d0)) 
            else  
              err = max(err,ABS(sum(i,j))) 
            endif  
          enddo  
          write(99,'(9999(1pE11.3))') (sum(i,j),j=1,N+1)
        enddo  
        write(99,*) "maximum error=",err
        print*,"matrix inversion error=",err
        stop
      endif  

      end
