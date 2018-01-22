************************************************************************
      subroutine INIT_TIMESTEP
************************************************************************
      use GRID,ONLY: N=>Npoints,zz,dt_diff_ex
      use STRUCT,ONLY: nHtot,Diff
      implicit none
      real*8 :: D,dt
      integer :: i

      dt = 9.D+99
      do i=2,N
        D  = 0.5*(Diff(i-1)+Diff(i)) 
        dt = MIN(dt,0.5*(zz(i)-zz(i-1))**2/D)
      enddo
      dt_diff_ex = dt
      print*,"explicit diffusion timestep=",dt_diff_ex

      end
