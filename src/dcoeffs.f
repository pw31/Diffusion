************************************************************************
      subroutine DCOEFFS
************************************************************************
      use GRID,ONLY: N=>Npoints,zz,d1l,d1m,d1r,d2l,d2m,d2r
      use PARAMETERS,ONLY: Hp
      use NATURE,ONLY: pi
      implicit none
      integer :: i
      real :: k,df,df2,h1,h2,x
      real,allocatable,dimension(:) :: hl,hr,f0,f1,f2
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
        k = 2.0*pi/Hp
        do i=1,N
          f0(i) = sin(k*zz(i))        ! test function
          f1(i) = cos(k*zz(i))*k        
          f2(i) =-sin(k*zz(i))*k**2
          !x = zz(i)/Hp 
          !f0(i) =  0.5*x**4 -1.0*x**3 +2.0*x**2 -3.0*x + 4.0
          !f1(i) = (2.0*x**3 -3.0*x**2 +4.0*x    -3.0)/Hp
          !f2(i) = (6.0*x**2 -6.0*x    +4.0          )/Hp**2
        enddo
        do i=1,N
          if (i==1) then
            df  = f0(1)*d1l(1) + f0(2)*d1m(1) + f0(3)*d1r(1) 
            df2 = 0.0
          else if (i==N) then
            df  = f0(N-2)*d1l(N) + f0(N-1)*d1m(N) + f0(N)*d1r(N) 
            df2 = 0.0
          else
            df  = f0(i-1)*d1l(i) + f0(i)*d1m(i) + f0(i+1)*d1r(i)
            df2 = f0(i-1)*d2l(i) + f0(i)*d2m(i) + f0(i+1)*d2r(i) 
          endif  
          print'(I4,99(1pE12.4))',i,f0(i),zz(i),
     >          f1(i),df,f1(i)-df,f2(i),df2,f2(i)-df2
        enddo
        stop
      endif  

      end
