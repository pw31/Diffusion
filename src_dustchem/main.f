      program DiffuDrift
      use PARAMETERS,ONLY: init,tsim,verbose
      use GRID,ONLY: dt_diff_ex
      use EXCHANGE,ONLY: chemcall,chemiter
      implicit none
      real*8  :: time,dt
      integer :: it,nout,next

      call INIT_NATURE
      call READ_PARAMETER
      call INIT_ELEMENTS1
      call INIT_CHEMISTRY
      call INIT_DUSTCHEM
      call READ_STRUCTURE
      call INIT_ELEMENTS2
      call INIT_GRID
      call INIT_DIFFUSION
      if (verbose>1) call OMP_SET_NUM_THREADS(1)

      nout = 0
      time = 0.d0
      dt   = 1.E-3*dt_diff_ex
      call INITIAL_CONDITIONS(nout,time,dt)
      next = nout
      if (next>=100)   next=next+1
      if (next>=300)   next=next+3
      if (next>=1000)  next=next+5
      if (next>=3000)  next=next+10
      if (next>=10000) next=next+30
      do 
        print* 
        print'("new timestep",i8,"  Dt=",1pE10.3," ...")',nout,dt 
        call DIFFUSION(time,0.5*dt,verbose)
        call SETTLING(time,dt,verbose)   ! changes dt
        call DIFFUSION(time,0.5*dt,verbose)
        call DUSTFORM(time,dt,verbose)
        time = time+dt
        nout = nout + 1
        if (nout>next) then
          call OUTPUT(nout,time,dt)
          next = nout
          if (next>=100)   next=next+1
          if (next>=300)   next=next+3
          if (next>=1000)  next=next+5
          if (next>=3000)  next=next+10
          if (next>=10000) next=next+30
        endif  
        if (time>tsim) exit
      enddo  

      print*
      print'("         smchem calls = ",I8)',chemcall
      print'("      iterations/call = ",0pF8.2)',
     >                     REAL(chemiter)/REAL(chemcall)

      end

      
