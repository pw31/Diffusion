************************************************************************
      SUBROUTINE  DUST_TSTEP(NN,yy,FF,tt,dt0,deltat,verbose,evap)
************************************************************************
      use NATURE,ONLY: pi,mic
      use PARAMETERS,ONLY: evap_model
      use EXCHANGE,ONLY: Fcall,Jcall,Jst,chi,ipoint
      use DUST_DATA,ONLY: NDUST
      use ELEMENTS,ONLY: NEPS,eps0,elnr,elnam
      implicit none
      integer,intent(in)  :: NN,verbose
      real*8,intent(inout):: yy(NN),FF(NN),tt
      real*8,intent(in)   :: dt0,deltat
      logical,intent(out) :: evap
      real*8  :: yold(NN)
      real*8  :: dt,t0,t1,t2,delt,dtold,rstep,bmix
      real*8  :: amean0,amean,tol=1.d-4
      integer :: i,it,ierr,nstep,nFcall,nJcall,nJac
      logical :: IS_NAN

* ----------------------------------------
* ***  variables for the LIMEX Solver  ***
* ----------------------------------------
      external :: FCN,JACOBI
      real*8   :: ropt(5)
      real*8   :: rtol(NN),atol(NN)
      integer  :: ipos(NN)
      integer  :: iopt(30),ifail(3)

      !----------------------------------------------
      ! ***  set parameters for the LIMEX-solver  ***
      !----------------------------------------------
      iopt(1)  = 0       ! how much output? (0=no, 1=standard, 2=more)
      iopt(2)  = 0       ! unit number for output (0:default=6)
      iopt(3)  = 0       ! solution output? (0=no)
      iopt(4)  = 0       ! unit number for solution output (0:default=6)
      iopt(5)  = 1       ! nonsigular matrix BB (0=singular, 1=nonsingular)
      iopt(6)  = 0       ! determination of initial values for FF,BB (1=yes)
      iopt(7)  = 0       ! analytic Jacobian? (0=numerical, 1=analytic)
      iopt(8)  = NN      ! Lower bandwidth of the Jacobian (Ndim=full)
      iopt(9)  = NN      ! Upper bandwidth of the Jacobian (Ndim=full)
      iopt(10) = 1       ! reuse of the Jacobian? (1=yes) 
      iopt(11) = 1       ! Switch for error toleranz (0=scalars, 1=vectors)
      iopt(12) = 0       ! Return after one integration time step? (0=no)
      iopt(13) = 0       ! Dense output option (0=no dense output)
      iopt(14) = 0       ! The number of equidistant dense output points
      iopt(15) = 0       ! unit for dense output
      iopt(16) = 0       ! type of call (0=first call, 1=continuation call)
      iopt(17) = 0       ! bevahior at t1 (0=exactly up to t1)
      iopt(18) = 0       ! Generation of of a PostScript plot of the Jacobian
      ropt(1)  = deltat  ! maximum allowed stepsize 
      ropt(2)  = 0.0     ! maximal distance of dense outputs (if iopt(13)=3)
      ropt(3)  = 0.0     ! upper limit for t (if iopt(17)=1)
      atol(:) = 1.d-299  ! absolute tolerance
      rtol(:) = tol      ! relative tolerance
      if (evap_model==1) then
        ipos(:) = 1      ! prevent YY(i)<0 if ipos(i)=1
      else   
        ipos(:) = 0
        ipos(1:4) = 1
        ipos(5+NDUST:NN) = 1
        atol(5:4+NDUST) = 1.d-40
      endif
  
      Fcall = 0
      Jcall = 0
      dt = MIN(dt0,1.E-3*deltat)
      delt = dt*30.d0
      if (yy(1)>0.d0) then
        amean0 = MIN((3.d0/(4.d0*pi))**(1.d0/3.d0)*yy(2)/yy(1),
     >               (3.d0/(4.d0*pi)*yy(4)/yy(1))**(1.d0/3.d0))
      else
        amean0 = 0.d0
      endif  
      evap = .false.

      do it=1,999999

        t0 = tt 
        t1 = tt
        t2 = min(deltat,tt+delt)
        yold = yy
        dtold = dt
#ifdef CHECK_NAN
        do i=1,NN
          if ((yy(i)<0.d0.and.ipos(i)==1).or.IS_NAN(yy(i))) then
            print*,"*** before limex",ipoint,yy 
            stop
          endif
        enddo  
#endif

        call LIMEX ( NN, FCN, JACOBI, t1, t2, yy, FF, 
     &               rtol, atol, dt, iopt, ropt, ipos, ifail )
        
#ifdef CHECK_NAN
        do i=1,NN
          if ((yy(i)<0.d0.and.ipos(i)==1).or.IS_NAN(yy(i))) then
            print*,"*** after limex",ipoint,t0,t1,t2
            print*,yold
            print*,yy 
            stop
          endif
        enddo  
#endif
        ierr   = ifail(1)
        nFcall = Fcall     !iopt(24)
        nJcall = iopt(25)
        nJac   = Jcall/NN  !iopt(29)
        nstep  = iopt(28)
        rstep  = (t1-t0)/(t2-t0)
        if (verbose>1) then
          print'(" ifail=",I3,"  calls=",I7,I6,"  rstep=",0pF8.6,
     >           "  t1/deltat=",1pE10.4,"  dt=",1pE8.2)',
     >           ierr,nFcall,nJac,rstep,t1/deltat,dtold
        endif  

        if (ierr==-46) then      ! too many step-size reductions
          tt = t1 
          delt = delt/10.d0
          dt = dtold/10.d0
          iopt(10) = 0           ! no re-use Jacobian
          iopt(16) = 0           ! new init call
          call DUST_RHS(NN,tt,yy,FF,.false.,ierr)
          evap = .true.
          do i=1,NDUST
            bmix = yy(4+i)/yy(4) 
            evap = evap.and.((chi(i)<0.d0).or.(bmix<1.e-3))
          enddo
          amean = MIN((3.d0/(4.d0*pi))**(1.d0/3.d0)*yy(2)/yy(1),
     >                (3.d0/(4.d0*pi)*yy(4)/yy(1))**(1.d0/3.d0))
          evap = evap.or.(amean<0.01*amean0) 
          if (evap) exit
        else if (ierr==-48) then ! too many integration steps
          tt = t1 
          dt = dtold
          iopt(10) = 1           ! re-use Jacobian
          iopt(16) = 0           ! new init call
        else if (ierr==-50) then ! step-size too small
          tt = t0
          yy = yold 
          delt = delt/10.d0
          dt = dtold/10.d0
          iopt(10) = 0           ! no re-use Jacobian
          iopt(16) = 0           ! new init call
        else if (ierr<0) then
          print*,"*** ierr=",ierr 
          stop 
        else  
          tt = t1
          delt = MIN(ropt(1),30.d0*dt)
          iopt(10) = 1           ! re-use Jacobian
          iopt(16) = 1           ! continuation call
        endif  

        if (dt<1.E-15*tt) then
          call DUST_RHS(NN,tt,yy,FF,.false.,ierr)
          evap = .true.
          do i=1,NDUST
            bmix = yy(4+i)/yy(4) 
            evap = evap.and.((chi(i)<0.d0).or.(bmix<1.e-3))
          enddo
          amean = MIN((3.d0/(4.d0*pi))**(1.d0/3.d0)*yy(2)/yy(1),
     >                (3.d0/(4.d0*pi)*yy(4)/yy(1))**(1.d0/3.d0))
          evap = evap.or.(amean<0.01*amean0) 
          if (evap) exit
          print*
          print*,"... tt+dt==tt problem."
          call DUST_RHS(NN,tt,yy,FF,.true.,ierr)
          stop
        endif
        if (tt.ge.deltat) exit

      enddo  

      end

