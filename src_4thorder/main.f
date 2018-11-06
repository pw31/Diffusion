      program DiffuDrift
      use PARAMETERS,ONLY: implicit,tindep,simple
      implicit none

      call INIT_NATURE
      call READ_PARAMETER
      if (simple) then
        call INIT_GRID
      else  
        call READ_STRUCTURE
        call FILLIN_GRID
      endif  
      call DCOEFFS
      call INITIAL_CONDITIONS
      if (tindep) then
        call DIFFUSION_TINDEPENDENT
      else if (implicit) then
        call DIFFUSION_IMPLICIT
      else
        call DIFFUSION_EXPLICIT
      endif  


      end

      
