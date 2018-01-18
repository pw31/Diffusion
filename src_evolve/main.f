
      program DiffuDrift

      use PARAMETERS,ONLY: init,tsim
      use GRID,ONLY: dt_diff_ex
      use EXCHANGE,ONLY: chemcall,chemiter,ieqcond,ieqconditer,
     >                   itransform
      use DATABASE,ONLY: NLAST
      implicit none
      real*8  :: time,dt,Tg,nHges,p,delta_z
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
      call INIT_CRUST

      nout = 0
      time = 0.d0
      dt   = 1.E-3*dt_diff_ex
      call INITIAL_CONDITIONS(nout,time,dt)
      next = nout

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
      
      !call DIFFUSION(time,dt,verbose)

      end

      
