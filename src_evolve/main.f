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

      !call DIFFUSION(time,0.5*dt,verbose)

      end

      
