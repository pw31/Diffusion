************************************************************************
      subroutine DCOEFFS
************************************************************************
      use GRID,ONLY: N=>Npoints,zz,d1l,d1m,d1r,d2l,d2m,d2r
      use PARAMETERS,ONLY: Hp,hmin,hmax
      implicit none
      integer :: i
      real*8 :: k,df,h1,h2
      real*8,allocatable,dimension(:) :: hl,hr,f0,f1,f2
      logical :: test=.false.

      allocate(hl(2:N),hr(1:N-1))
      allocate(d1l(1:N),d1m(1:N),d1r(1:N))
      allocate(d2l(1:N),d2m(1:N),d2r(1:N))
      allocate(f0(N),f1(N),f2(N))

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
        k = 3.0/(Hp*(hmax-hmin))
        do i=1,N
          f0(i) = sin(k*zz(i))        ! test function
          f1(i) = cos(k*zz(i))*k        
          f2(i) =-sin(k*zz(i))*k**2
          f0(i) = 2.0*zz(i)**2-zz(i)+0.5
          f1(i) = 4.0*zz(i)-1.0
          f2(i) = 4.0
        enddo
        do i=2,N-1
          df = f0(i-1)*d1l(i) + f0(i)*d1m(i) + f0(i+1)*d1r(i)
          print'(I4,2(0pF12.6),0pF10.6)',i,f1(i),df,f1(i)-df
        enddo
        do i=2,N-1
          df = f0(i-1)*d2l(i) + f0(i)*d2m(i) + f0(i+1)*d2r(i)
          print'(I4,2(0pF12.6),0pF10.6)',i,f2(i),df,f2(i)-df
        enddo
        df = f0(1)*d1l(1) + f0(2)*d1m(1) + f0(3)*d1r(1) 
        print'(I4,2(0pF12.6),0pF10.6)',1,f1(1),df,f1(1)-df
        df = f0(N-2)*d1l(N) + f0(N-1)*d1m(N) + f0(N)*d1r(N) 
        print'(I4,2(0pF12.6),0pF10.6)',N,f1(N),df,f1(N)-df
        stop
      endif  

      end
