
      program DiffuDrift

      use NATURE,ONLY: km
      use PARAMETERS,ONLY: init,tsim,verbose
      use GRID,ONLY: dt_diff_ex,zz
      use STRUCT,ONLY: Temp,nHtot,press,nHeps,crust_depth,crust_beta
      use ELEMENTS,ONLY: NELEM,NEPS,eps0,elnr
      use CHEMISTRY,ONLY: NMOLE
      use DUST_DATA,ONLY: NDUST,dust_vol,dust_nam
      use EXCHANGE,ONLY: chemcall,chemiter,ieqcond,ieqconditer,
     >                   itransform,inactive,nmol,chi
      use DATABASE,ONLY: NLAST
      implicit none
      real*8  :: time,dt,Tg,nHges,p,delta_z
      integer :: i,j,el,it,nout,next
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp),allocatable,dimension(:) :: eps,Sat,eldust
      real(kind=qp) :: sum_beta

      call INIT_NATURE
      call READ_PARAMETER
      call INIT_ELEMENTS1
      call INIT_CHEMISTRY
      call INIT_DUSTCHEM
      call READ_STRUCTURE
      call INIT_ELEMENTS2
      call INIT_GRID
      call INIT_DIFFUSION

      allocate(eps(NELEM),Sat(NDUST),eldust(NDUST))
      allocate(nmol(NMOLE),chi(NDUST),inactive(NMOLE))
      i = 1
      Tg      = Temp(i)
      p       = press(i) 
      nHges   = nHtot(i)
      inactive = .false.
      call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)

      delta_z = zz(2)-zz(1)
      crust_depth = 1.d0*km
      crust_beta(:) = 0.Q0
      sum_beta = 0.Q0
      do i=1,NDUST
        if (eldust(i).le.0.Q0) cycle
        crust_beta(i) = eldust(i)*dust_vol(i)
        sum_beta = sum_beta + crust_beta(i)
      enddo  
      crust_beta = crust_beta/sum_beta
      do i=1,NDUST
        if (eldust(i).le.0.Q0) cycle
        print*,dust_nam(i),REAL(crust_beta(i))
      enddo  
      print*,delta_z,crust_depth

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

      
