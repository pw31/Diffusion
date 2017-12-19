      program DiffuDrift
      use PARAMETERS,ONLY: init,tsim,verbose
      use GRID,ONLY: dt_diff_ex
      implicit none
      real*8  :: time,dt
      integer :: it,nout

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
      do 
        print* 
        print'("new timestep",i8,"  Dt=",1pE10.3," ...")',nout,dt 
        call DIFFUSION(time,0.5*dt,verbose)
        call SETTLING(time,dt,verbose)   ! changes dt
        call DIFFUSION(time,0.5*dt,verbose)
        call DUSTFORM(time,dt,verbose)
        time = time+dt
        call OUTPUT(nout,time,dt)
        if (time>tsim) exit
      enddo  

      end

      
