************************************************************************
      subroutine TIMESTEP(dt_settle,dt_dustform,dt)
************************************************************************
      use PARAMETERS,ONLY: dtfix,precision
      implicit none
      real*8,intent(in) :: dt_settle,dt_dustform
      real*8,intent(inout) :: dt

      !print*,"TIMESTEP:",dt_settle,dt_dustform,dt
      if (dtfix>0.0) then
        dt = dtfix
      else  
        if (MIN(dt_settle,dt_dustform)>1.4*dt) then
          dt = dt*1.3
          print'("... increasing timestep =",1pE9.2,
     >           " settling =",1pE9.2," dust formation =",1pE9.2)',
     >         dt,dt_settle,dt_dustform
        else if (MIN(dt_settle,dt_dustform)<dt) then
          dt = dt*0.7
          print'("... decreasing timestep =",1pE9.2,
     >           " settling =",1pE9.2," dust formation =",1pE9.2)',
     >         dt,dt_settle,dt_dustform
        endif
      endif

      end
          
          
