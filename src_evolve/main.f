
      program DiffuDrift

      use PARAMETERS,ONLY: tsim,dt_init,dt_increase,dt_max,verbose
      use GRID,ONLY: dt_diff_ex
      use EXCHANGE,ONLY: chemcall,chemiter,ieqcond,ieqconditer,
     >                   itransform
      use DATABASE,ONLY: NLAST
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
      call INIT_TIMESTEP
      call INIT_CRUST

      nout = 0
      time = 0.d0
      dt   = dt_init
      call INITIAL_CONDITIONS(nout,time,dt)
      call OUTPUT(nout,time,dt)
      next = nout

      do it=1,999999
        print* 
        call WARM_UP(time)
        print'("new timestep",i8,"  Dt=",1pE10.3," ...")',nout,dt 
        call DIFFUSION(time,dt,verbose)
        call UPDATE_CRUST
        time = time+dt
        nout = nout + 1
        if (nout>next) then
          call OUTPUT(nout,time,dt)
          next = nout
          dt = MIN(dt_max,dt*dt_increase)
        endif  
        if (time>tsim) exit
      enddo  
      call OUTPUT(nout,time,dt)

      print*
      print'("         smchem calls = ",I8)',chemcall
      print'("      iterations/call = ",0pF8.2)',
     >                     REAL(chemiter)/REAL(chemcall)
      print'("eq condensation calls = ",I8)',ieqcond
      print'("   eq iterations/call = ",0pF8.2)',
     >                  REAL(ieqconditer)/REAL(ieqcond)
      print'("      transform calls = ",I8)',itransform
      NLAST=0         ! also save replaced database entries
      call SAVE_DBASE
      
      end

      
