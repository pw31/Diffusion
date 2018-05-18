**********************************************************************
      SUBROUTINE WARM_UP(time)
**********************************************************************
      use PARAMETERS,ONLY: heatrate,verbose
      use GRID,ONLY: Npoints
      use STRUCT,ONLY: Temp,Temp0
      implicit none
      real*8,intent(in) :: time
      integer :: i

      do i=-2,Npoints
        Temp(i) = Temp0(i) + heatrate*time
      enddo
      if (verbose<=0) return
      print*
      print'(" SURFACE TEMPERATURE =",0pF8.2," K")',Temp(0)
      print'(A32)'," ==============================="

      end
