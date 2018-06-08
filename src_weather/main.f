      program DiffuDrift
      use PARAMETERS,ONLY: init,tsim,verbose,dtfix,dnfix
      use GRID,ONLY: dt_diff_ex
      use EXCHANGE,ONLY: chemcall,chemiter
      implicit none
      real*8  :: time,dt,dt_settle,dt_dustform
      integer :: it,nout,next

      call INIT_NATURE
      call READ_PARAMETER
      call INIT_ELEMENTS1
      call INIT_CHEMISTRY
      call INIT_DUSTCHEM
      !call READ_STRUCTURE_PHOENIX
      call READ_STRUCTURE_VPARM(-135,0)   ! night side
      call INIT_ELEMENTS2
      call INIT_GRID
      call INIT_DIFFUSION
      if (verbose>1) call OMP_SET_NUM_THREADS(1)

      nout = 0
      time = 0.d0
      dt   = 1.E-3*dt_diff_ex
      call INITIAL_CONDITIONS(nout,time,dt)
      if (dtfix>0.0) dt=dtfix
      next = nout
      if (next>=100) next=next+1
      if (next>=300) next=next+3
      if (next>=1000) next=next+5
      if (next>=3000) next=next+10
      if (next>=10000) next=next+30
      if (next>=30000) next=next+50
      if (dnfix>0) next=nout-1+dnfix
      do it=1,9999999
        print* 
        print'("new timestep",i8,"  Dt=",1pE10.3," ...")',nout,dt 
        call DIFFUSION(time,0.5*dt,verbose)
        call SETTLING (time,1.0*dt,dt_settle,verbose)   
        call DIFFUSION(time,0.5*dt,verbose)
        call DUSTFORM (time,1.0*dt,dt_dustform,verbose)
        time = time+dt
        if (time>tsim) exit
        nout = nout + 1
        if (nout>next) then
          call OUTPUT(nout,time,dt)
          next = nout
          if (next>=100) next=next+1
          if (next>=300) next=next+3
          if (next>=1000) next=next+5
          if (next>=3000) next=next+10
          if (next>=10000) next=next+30
          if (next>=30000) next=next+50
          if (dnfix>0) next=nout-1+dnfix
        endif  
        call TIMESTEP(dt_settle,dt_dustform,dt)
      enddo  

      print*
      print'("         smchem calls = ",I12)',chemcall
      print'("      iterations/call = ",0pF12.3)',
     >                     REAL(chemiter)/REAL(chemcall)

      end

      
