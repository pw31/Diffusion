************************************************************************
      subroutine DCOEFFS
************************************************************************
      use GRID,ONLY: N=>Npoints,zz,d1l2,d1l1,d1m,d1r1,d1r2,d2l2,d2l1,
     >d2m,d2r1,d2r2
      use PARAMETERS,ONLY: Hp
      implicit none
      integer :: i
      real*8 :: k,x,h,df,df2
      real*8,allocatable,dimension(:) :: hl1, hl2, hr1, hr2,f0,f1,f2
      logical :: test=.true.

      allocate(hl2(1:N),hl1(1:N),hr1(1:N),hr2(1:N))
      allocate(d1l2(1:N),d1l1(1:N),d1m(1:N),d1r1(1:N),d1r2(1:N))
      allocate(d2l2(1:N),d2l1(1:N),d2m(1:N),d2r1(1:N),d2r2(1:N))
      allocate(f0(N),f1(N),f2(N))

      !--- compute grid point differences ---
      do i=3,N
        hl1(i) = zz(i)-zz(i-1)
        hl2(i) = zz(i)-zz(i-2)
      enddo
      do i=1,N-2
        hr1(i) = zz(i+1)-zz(i)
        hr2(i) = zz(i+2)-zz(i)
      enddo

      !--- compute 1st and 2nd derivative coefficients ---
      do i=3,N-2 

        !d1l(i) = -hr(i)/((hr(i)+hl(i))*hl(i))
        !d1m(i) =  (hr(i)-hl(i))/(hl(i)*hr(i))
        !d1r(i) =  hl(i)/((hr(i)+hl(i))*hr(i))
        d1l2(i) = -(hl1(i)*hr1(i)*hr2(i))/((hl1(i)-hl2(i))*hl2(i)*
     >(hl2(i)+hr1(i))*(hl2(i)+hr2(i)))
        d1l1(i) = (hl2(i)*hr1(i)*hr2(i))/(hl1(i)*(hl1(i)-hl2(i))*
     >(hl1(i)+hr1(i))*(hl1(i)+hr2(i)))
        d1m(i)  = 1/(hl1(i))+1/(hr2(i))-(hr1(i)+hr2(i))/(hr1(i)*
     >hr2(i))
        d1r1(i) = (hl1(i)*hl2(i)*hr2(i))/(hr1(i)*(hl1(i)+hr1(i))*
     >(hl2(i)+hr1(i))*(hr2(i)-hr1(i)))
        d1r2(i) = (hl1(i)*hl2(i)*hr1(i))/((hr1(i)-hr2(i))*hr2(i)*
     >(hl1(i)+hr2(i))*(hl2(i)+hr2(i)))

        !d2l(i) =  2.0/((hr(i)+hl(i))*hl(i))
        !d2m(i) = - 2.0/(hr(i)*hl(i))
        !d2r(i) =  2.0/((hr(i)+hl(i))*hr(i))
        d2l2(i) = (-2*hr1(i)*hr2(i)+2*hl1(i)*(hr1(i)+hr2(i)))/
     >((hl1(i)-hl2(i))*hl2(i)*(hl2(i)+hr1(i))*(hl2(i)+hr2(i)))
        d2l1(i) = (2*hr1(i)*hr2(i)-2*hl2(i)*(hr1(i)+hr2(i)))/(hl1(i)*
     >(hl1(i)-hl2(i))*(hl1(i)+hr1(i))*(hl1(i)+hr2(i)))
        d2m(i)  = (2*(hl1(i)*(hl2(i)-hr1(i)-hr2(i))+hr1(i)*hr2(i)-
     >hl2(i)*(hr1(i)+hr2(i))))/(hl1(i)*hl2(i)*hr1(i)*hr2(i))
        d2r1(i) = (2*hl1(i)*hl2(i)-2*(hl1(i)+hl2(i))*hr2(i))/(hr1(i)*
     >(hl1(i)+hr1(i))*(hl2(i)+hr1(i))*(hr1(i)-hr2(i)))
        d2r2(i) = (-2*hl1(i)*hl2(i)+2*(hl1(i)+hl2(i))*hr1(i))/((hr1(i)-
     >hr2(i))*hr2(i)*(hl1(i)+hr2(i))*(hl2(i)+hr2(i)))

      enddo

      !--- compute 1st and 2nd derivative coefficients at inner point of boundaries ---
      !d1m(1) = -1.d0/hr(1)
      !d1r(1) =  1.d0/hr(1)  ! first order
      hl2(2)   = zz(1) - zz(2)
      hl1(2)   = zz(3) - zz(2)
      hr1(2)   = zz(4) - zz(2)
      hr2(2)   = zz(5) - zz(2)

      d1l2(2)  =(hl1(2)*hr1(2)*hr2(2))/((hl1(2) - hl2(2))*hl2(2)*
     >(-hr1(2) + hl2(2))*(-hr2(2) + hl2(2)))
      d1l1(2)  =(-1.0)/(hl1(2)) -(1.0)/(hr1(2)) - (hr2(2)+hl2(2))/
     >(hr2(2)*hl2(2))
      d1m(2)   =(hr1(2)*hr2(2)*hl2(2))/(hl1(2)*(-hl1(2) + hr1(2))*
     >(-hl1(2) + hr2(2))*(-hl1(2) + hl2(2)))
      d1r1(2)  =(hl1(2)*hr2(2)*hl2(2))/((hl1(2) - hr1(2))*hr1(2)*
     >(hr1(2) - hr2(2))*(hr1(2) - hl2(2)))
      d1r2(2)  =(hl1(2)*hr1(2)*hl2(2))/((hl1(2) - hr2(2))*hr2(2)*
     >(-hr1(2) + hr2(2))*(hr2(2) - hl2(2)))

      d2l2(2)  =-(2.0*(hr1(2)*hr2(2) + hl1(2)*(hr1(2) + hr2(2))))/
     >((hl1(2) - hl2(2))*hl2(2)*(-hr1(2) + hl2(2))*(-hr2(2) + hl2(2)))
      d2l1(2)  =(2.0*(hr2(2)*hl2(2) + hr1(2)*(hr2(2) + hl2(2)) +
     >hl1(2)*(hr1(2) + hr2(2) + hl2(2))))/(hl1(2)*hr1(2)*hr2(2)*hl2(2))
      d2m(2)   =(2.0*(hr2(2)*hl2(2) + hr1(2)*(hr2(2) + hl2(2))))/
     >(hl1(2)*(hl1(2) - hr1(2))*(hl1(2) - hr2(2))*(hl1(2) - hl2(2)))
      d2r1(2)  =-(2.0*(hr2(2)*hl2(2) + hl1(2)*(hr2(2) + hl2(2))))/
     >((hl1(2) - hr1(2))*hr1(2)*(hr1(2) - hr2(2))*(hr1(2) - hl2(2)))
      d2r2(2)  =-(2.0*(hr1(2)*hl2(2) + hl1(2)*(hr1(2) + hl2(2))))/
     >((hl1(2) - hr2(2))*hr2(2)*(-hr1(2) + hr2(2))* (hr2(2) - hl2(2)))

      !d1l(N) = -1.d0/hl(N)
      !d1m(N) =  1.d0/hl(N)
      hl2(N-1)   = zz(N-4) - zz(N-1)
      hl1(N-1)   = zz(N-3) - zz(N-1)
      hr1(N-1)   = zz(N-2) - zz(N-1)
      hr2(N-1)   = zz(N) - zz(N-1)
      d1l2(N-1)  = (hr1(N-1)*hl1(N-1)*hr2(N-1))/((hr1(N-1) -
     >hl2(N-1))*hl2(N-1)*
     >(-hl1(N-1) + hl2(N-1))*(hl2(N-1) - hr2(N-1)))
      d1l1(N-1)  =(hr1(N-1)*hl2(N-1)*hr2(N-1))/((hr1(N-1)-hl1(N-1))*
     >hl1(N-1)*(hl1(N-1) - hl2(N-1))*(hl1(N-1) - hr2(N-1)))
      d1m(N-1)   =(hl1(N-1)*hl2(N-1)*hr2(N-1))/(hr1(N-1)*(-hr1(N-1) +
     >hl1(N-1))*(-hr1(N-1) + hl2(N-1))*(-hr1(N-1) + hr2(N-1)))
      d1r1(N-1)  =-(1.0)/(hr1(N-1))-(1.0)/(hl1(N-1))-(hl2(N-1)+
     >hr2(N-1))/(hl2(N-1)*hr2(N-1))
      d1r2(N-1)  =(hr1(N-1)*hl1(N-1)*hl2(N-1))/((hr1(N-1) - hr2(N-1))*
     >hr2(N-1)*(-hl1(N-1) + hr2(N-1))*(-hl2(N-1) + hr2(N-1)))

      d2l2(N-1)  =-(2.0*(hl1(N-1)*hr2(N-1) + hr1(N-1)*(hl1(N-1) +
     >hr2(N-1))))/((hr1(N-1) - hl2(N-1))*hl2(N-1)*(-hl1(N-1) +
     >hl2(N-1))*(hl2(N-1) - hr2(N-1)))
      d2l1(N-1)  =-(2.0*(hl2(N-1)*hr2(N-1) + hr1(N-1)*(hl2(N-1) +
     >hr2(N-1))))/((hr1(N-1) - hl1(N-1))*hl1(N-1)*(hl1(N-1) -
     >hl2(N-1))*(hl1(N-1) - hr2(N-1)))
      d2m(N-1)   =(2.0*(hl2(N-1)*hr2(N-1) + hl1(N-1)*(hl2(N-1) +
     >hr2(N-1))))/(hr1(N-1)*(hr1(N-1) - hl1(N-1))*(hr1(N-1) -
     >hl2(N-1))*(hr1(N-1) - hr2(N-1)))
      d2r1(N-1)  =(2.0*(hl2(N-1)*hr2(N-1) + hl1(N-1)*(hl2(N-1) +
     >hr2(N-1)) + hr1(N-1)*(hl1(N-1)+hl2(N-1)+hr2(N-1))))/(
     >hr1(N-1)*hl1(N-1)*hl2(N-1)*hr2(N-1))
      d2r2(N-1)  =-(2.0* (hl1(N-1)*hl2(N-1) + hr1(N-1)*(hl1(N-1) +
     >hl2(N-1))))/((hr1(N-1) - hr2(N-1))*hr2(N-1)*(-hl1(N-1) +
     >hr2(N-1))*(-hl2(N-1) + hr2(N-1)))


      !--- compute 1st derivative coefficients at very boundaries ---
      !d1m(1) = -1.d0/hr(1)
      !d1r(1) =  1.d0/hr(1)  ! first order
      hl2(1)   = zz(2) - zz(1)
      hl1(1)   = zz(3) - zz(1)
      hr1(1)   = zz(4) - zz(1)
      hr2(1)   = zz(5) - zz(1)
      d1l2(1)  = -(1.0)/(hl2(1))-(1.0)/(hl1(1))-(hr1(1)+hr2(1))/
     >(hr1(1)*hr2(1))
      d1l1(1)  = (hl1(1)*hr1(1)*hr2(1))/(hl2(1)*(-hl2(1) + hl1(1))*
     >(-hl2(1) + hr1(1))*(-hl2(1) + hr2(1)))
      d1m(1)   = (hl2(1)*hr1(1)*hr2(1))/((hl2(1) - hl1(1))*hl1(1)*
     >(hl1(1) - hr1(1))*(hl1(1) - hr2(1)))
      d1r1(1)  = (hl2(1)*hl1(1)*hr2(1))/((hl2(1) - hr1(1))*hr1(1)*
     >(-hl1(1) + hr1(1))*(hr1(1) - hr2(1)))
      d1r2(1)  = (hl2(1)*hl1(1)*hr1(1))/((hl2(1) - hr2(1))*hr2(1)*
     >(-hl1(1) + hr2(1))*(-hr1(1) + hr2(1)))
      !d1l(N) = -1.d0/hl(N)
      !d1m(N) =  1.d0/hl(N)
      hl2(N)   = zz(N-4) - zz(N)
      hl1(N)   = zz(N-3) - zz(N)
      hr1(N)   = zz(N-2) - zz(N)
      hr2(N)   = zz(N-1) - zz(N)
      d1l2(N)  = (hr2(N)*hr1(N)*hl1(N))/((hr2(N) - hl2(N))*hl2(N)*
     >(-hr1(N) + hl2(N))*(-hl1(N) + hl2(N)))
      d1l1(N)  = (hr2(N)*hr1(N)*hl2(N))/((hr2(N) - hl1(N))*hl1(N)*
     >(-hr1(N) + hl1(N))*(hl1(N) - hl2(N)))
      d1m(N)   = (hr2(N)*hl1(N)*hl2(N))/((hr2(N) - hr1(N))*hr1(N)*
     >(hr1(N) - hl1(N))*(hr1(N) - hl2(N)))
      d1r1(N)  = (hr1(N)*hl1(N)*hl2(N))/(hr2(N)*(-hr2(N) + hr1(N))*
     >(-hr2(N) + hl1(N))*(-hr2(N) + hl2(N)))
      d1r2(N)  = -(1.0)/(hr2(N))-(1.0)/(hr1(N))-(hl1(N)+hl2(N))/
     >(hl1(N)*hl2(N))


      if (test) then
        !--- test derivatives ---
        k = 3.0/Hp
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
            df = f0(i)*d1l2(i) + f0(i+1)*d1l1(i) + f0(i+2)*d1m(i) +
     >           f0(i+3)*d1r1(i) + f0(i+4)*d1r2(i)
            df2 = 0.0
          else if (i==2) then
            df = f0(i-1)*d1l2(i) + f0(i)*d1l1(i) + f0(i+1)*d1m(i) +
     >           f0(i+2)*d1r1(i) + f0(i+3)*d1r2(i)
            df2= f0(i-1)*d2l2(i) + f0(i)*d2l1(i) + f0(i+1)*d2m(i) +
     >           f0(i+2)*d2r1(i) + f0(i+3)*d2r2(i)
          else if (i==N-1) then
            df = f0(i-3)*d1l2(i) + f0(i-2)*d1l1(i) + f0(i-1)*d1m(i) +
     >           f0(i)*d1r1(i) + f0(i+1)*d1r2(i)
            df2= f0(i-3)*d2l2(i) + f0(i-2)*d2l1(i) + f0(i-1)*d2m(i) +
     >           f0(i)*d2r1(i) + f0(i+1)*d2r2(i)
          else if (i==N) then
            df = f0(i-4)*d1l2(i) + f0(i-3)*d1l1(i) + f0(i-2)*d1m(i) +
     >           f0(i-1)*d1r1(i) + f0(i)*d1r2(i)
            df2= 0.0
          else
            df = f0(i-2)*d1l2(i) + f0(i-1)*d1l1(i) + f0(i)*d1m(i) +
     >           f0(i+1)*d1r1(i) + f0(i+2)*d1r2(i)
            df2= f0(i-2)*d2l2(i) + f0(i-1)*d2l1(i) + f0(i)*d2m(i) +
     >           f0(i+1)*d2r1(i) + f0(i+2)*d2r2(i)
          end if
          if (i==1) then
            h = zz(i+1)-zz(i) 
          else
            h = zz(i)-zz(i-1) 
          endif
          print'(I4,99(1pE13.4))',i,f0(i),zz(i),
     >          f1(i),df,f1(i)-df,f2(i),df2,f2(i)-df2
          !print'(I4,99(1pE20.12))',i,
     >    !      d1l2(i)*h,d1l1(i)*h,d1m(i)*h,d1r1(i)*h,d1r2(i)*h

        enddo
        !check if the calculation in the boundary are correct
        !print'(I4,2(0pF12.6),0pF10.6)',2,f1(1),df,f1(1)-df
        !f0(i-1)*d2l1(i) +x f0(i)*d2m(i) + f0(i+1)*d2r1(i)
        !print'(I4,2(0pF12.6),0pF10.6)',i,f2(i),df,f2(i)-df
        !enddo
        !df = f0(1)*d1l1(1) + f0(2)*d1m(1) + f0(3)*d1r1(1)
        !print'(I4,2(0pF12.6),0pF10.6)',1,f1(1),df,f1(1)-df
        !df = f0(N-2)*d1l1(N) + f0(N-1)*d1m(N) + f0(N)*d1r1(N)
        !print'(I4,2(0pF12.6),0pF10.6)',N,f1(N),df,f1(N)-df
        stop
      endif  

      end
