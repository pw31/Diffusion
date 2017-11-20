************************************************************************
      subroutine INITIAL_CONDITIONS
************************************************************************
      use NATURE,ONLY: pi
      use PARAMETERS,ONLY: tnull,init
      use GRID,ONLY: N=>Npoints,xx,zz
      use STRUCT,ONLY: nHtot,Diff
      implicit none
      integer :: i,nn
      real*8 :: L,k,ww,z0

      allocate(xx(N))
      tnull = 0.d0
      if (init==1) then
        xx(:) = 1.d0        ! constant concentration
      else if (init==2) then
        xx(:) = 0.d0
        xx(1) = 1.d0        ! inflow from below
      else if (init==3) then  
        !--- this one has a simple analytic solution ---
        nn = 1
        L  = zz(N)-zz(1)
        k  = pi/(nn*L)
        do i=1,N
          xx(i) = sin(k*(zz(i)-zz(1)))  
        enddo  
      else if (init==4) then 
        !--- this one has an analytic solution, too ---
        tnull = 0.0002
        ww = 2.d0*SQRT(nHtot(1)*Diff(1)*tnull)
        z0 = 0.5*(zz(N)+zz(1))
        do i=1,N
          xx(i) = exp(-((zz(i)-z0)/ww)**2)
        enddo  
      else if (init==5) then  
        xx(:) = 0.d0
        xx(N) = 1.d0        ! inflow from above
      else
        print*,"*** init=",init," not recognised."
        stop
      endif
      end
