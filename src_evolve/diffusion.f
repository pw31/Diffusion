************************************************************************
      subroutine DIFFUSION(time,deltat,verbose)
************************************************************************
      use PARAMETERS,ONLY: implicit
      use GRID,ONLY: zz,xlower,xupper,dt_diff_ex
      use STRUCT,ONLY: crust_Ncond,crust_gaseps,Temp,Diff,nHtot,nHeps      
      use ELEMENTS,ONLY: NELEM,elnam
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use EXCHANGE,ONLY: nat,nmol
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in) :: time,deltat
      integer,intent(in) :: verbose
      real(kind=qp) :: eps(NELEM),Sat(NDUST)
      real*8 :: nH,Tg,flux,fmin,dz,first_eps,second_eps
      logical :: evaporates(NDUST),limiting(NELEM)
      integer :: s,i,el,emin
      integer :: Nlim(NELEM),Dlim(NELEM,10)

      xlower = crust_gaseps

      !------------------------------------------------
      ! ***  Consider layer=2 to determine whether  ***
      ! ***  solids grow or evaporate.              ***
      ! ***  Determine flux limiting elements.      ***
      !------------------------------------------------
      Tg = Temp(1)
      nH = nHtot(1)
      eps(:) = nHeps(:,2)/nHtot(2)
      call GGCHEM(nH,Tg,eps,.false.,0)
      call SUPERSAT(Tg,nat,nmol,Sat)
      evaporates(:) = .false.
      limiting(:) = .false.
      Nlim(:) = 0
      do s=1,NDUST
        if (crust_Ncond(s)<=0.Q0) cycle
        if (Sat(s)<1.Q0) evaporates(s)=.true. 
        print'(A15,"  evaporates=",L1)',dust_nam(s),evaporates(s)
        dz = zz(2)-zz(1)
        fmin = 9.E+99
        emin = 0
        do i=1,dust_nel(s)
          el = dust_el(s,i)
          flux = -Diff(1)*nHtot(1)*(eps(el)-crust_gaseps(el))/dz
          if (ABS(flux)/dust_nu(s,i)<fmin) then
            fmin = ABS(flux)/dust_nu(s,i) 
            emin = el
          endif   
          print'(A3," flux=",1pE12.4," stoich=",I2)',
     >          elnam(el),flux,dust_nu(s,i)
        enddo
        print'("limiting element =",A3)',elnam(emin)
        limiting(emin) = .true.
        Nlim(emin) = Nlim(emin)+1
        Dlim(emin,Nlim(emin)) = s
      enddo

      if (implicit.and.deltat>100*dt_diff_ex) then
        if (verbose>1) then
          print*
          print*,"entering IMPLICIT DIFFUSION ..."
          print*,"==============================="
        endif  
        call DIFFUSION_IMPLICIT(time,deltat,limiting,Nlim,Dlim,verbose)
      else
        if (verbose>1) then
          print*
          print*,"entering EXPLICIT DIFFUSION ..."
          print*,"==============================="
        endif  
        call DIFFUSION_EXPLICIT(time,deltat,limiting,Nlim,Dlim,verbose)
      endif
      end


************************************************************************
      subroutine DIFFUSION_EXPLICIT(time0,deltat,limiting,Nlim,Dlim,
     >                              verbose)
************************************************************************
      use NATURE,ONLY: pi
      use GRID,ONLY: N=>Npoints,zz,d1l,d1m,d1r,d2l,d2m,d2r,dt_diff_ex,
     >               xlower,xupper
      use PARAMETERS,ONLY: dust_diffuse,bc_low,bc_high,
     >                     outflux,inrate,outrate,vin,vout
      use STRUCT,ONLY: Diff,nHtot,nHeps,crust_Neps,crust_Ncond
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use ELEMENTS,ONLY: NELEM,elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      implicit none
      real*8,intent(IN) :: time0,deltat
      logical,intent(in) :: limiting(NELEM)
      integer,intent(in) :: Nlim(NELEM),Dlim(NELEM,10),verbose
      integer :: i,it,e,el,round
      real*8,dimension(N) :: xx,rate
      real*8 :: influx(NELEM)
      real*8 :: D,nD,d1,d2,d1nD,time,dt
      character :: CR = CHAR(13)
      logical :: IS_NAN

      time = 0.d0
      do it=1,9999999

        dt = MIN(dt_diff_ex,deltat-time)

        do round=1,2

          if (round==2) then
            !---------------------------------------------------------
            ! ***  figure out influx of the non-limiting elements  ***
            ! ***  from the influxes of the limiting elements and  ***
            ! ***  the stoichiometry of the crust                  ***
            !---------------------------------------------------------
            call INFLUXES(influx,limiting,Nlim,Dlim) 
          endif
   
          do e=1,NELM
            if (e==iel) cycle
            el = elnum(e)
            if (round==1.and..not.limiting(el)) cycle
            if (round==2.and.limiting(el)) cycle

            xx(:) = nHeps(el,:)/nHtot(:) 
            xx(1) = xlower(el)            ! (not sure ...)

            !-------------------------------------------
            ! ***  d/dt(nH*x) = d/dz(nH*Diff*dx/dz)  ***
            !-------------------------------------------
            rate(:) = 0.0
            do i=2,N-1
              nD   = nHtot(i)*Diff(i)  
              d1   = d1l(i)*xx(i-1) + d1m(i)*xx(i) + d1r(i)*xx(i+1)
              d2   = d2l(i)*xx(i-1) + d2m(i)*xx(i) + d2r(i)*xx(i+1)
              d1nD = d1l(i)*nHtot(i-1)*Diff(i-1) 
     >             + d1m(i)*nHtot(i)  *Diff(i) 
     >             + d1r(i)*nHtot(i+1)*Diff(i+1) 
              rate(i) = nD*d2 + d1nD*d1
            enddo

            !----------------------------
            ! ***  explicit timestep  ***
            !----------------------------
            do i=2,N-1
              xx(i) = xx(i) + rate(i)*dt/nHtot(i)
            enddo  

            !------------------------------
            ! ***  boundary conditions  ***
            !------------------------------
            nD = nHtot(1)*Diff(1)  
            if (limiting(el)) then   
              xx(1) = xlower(el)              ! constant concentration
              influx(el) = -nD
     >              *(d1l(1)*xx(1) + d1m(1)*xx(2) + d1r(1)*xx(3))
              print'(A3,"  influx=",1pE12.4)',elnam(el),influx(el) 
            else                              
              xx(1) = (-influx(el)/nD - d1m(1)*xx(2) - d1r(1)*xx(3))
     >               /d1l(1)                  ! constant flux 
            endif  
            nD = nHtot(N)*Diff(N)  
            if (bc_high==1) then
              xx(N) = xupper(el)              ! const concentration
              outflux = -nD
     >                *(d1l(N)*xx(N-2) + d1m(N)*xx(N-1) + d1r(N)*xx(N)) 
            else if (bc_high==2) then   
              xx(N) = (-outflux/nD - d1l(N)*xx(N-2) - d1m(N)*xx(N-1))
     >                /d1r(N)                 ! constant flux 
            else if (bc_high==3) then   
              outflux = nHtot(N)*xx(N)*outrate*vout  
              xx(N) = -(d1l(N)*xx(N-2) + d1m(N)*xx(N-1))
     >                /(d1r(N) + outrate*vout/Diff(N))
            endif   

            !------------------------------------------
            ! ***  map solution on atmosphere grid  ***
            !------------------------------------------
            do i=1,N
              nHeps(el,i) = nHtot(i)*xx(i)
              if (xx(i)<0.d0) then
                print*,i
                print*,xx
                stop "should not occur." 
                if (i==N) then 
                  nHeps(el,i) = nHtot(i)*xx(N-1)
                else  
                  nHeps(el,i) = nHtot(i)*1.E-50
                endif
              endif    
            enddo 

          enddo
        enddo  

        do e=1,NELM
          if (e==iel) cycle
          el = elnum(e)
          !----------------------------------------
          ! ***  update crust column densities  ***
          !----------------------------------------
          if (crust_Neps(el)>0.Q0) then   
            print*,elnam(el),REAL(crust_Neps(el)),influx(el)*dt
            crust_Neps(el) = crust_Neps(el) - influx(el)*dt
            if (crust_Neps(el)<0.Q0) then
              print*,elnam(el),"negative crust column density"
              stop
            endif  
          endif  

        enddo  

        time = time + dt
        if (verbose>0) write(*,'(TL10," DIFFUSION:",I8,A,$)') it,CR
        if (time.ge.deltat) exit

      enddo  
      if (verbose>0) print'(" DIFFUSION:",I8,"  time=",1pE11.4,
     >               "  Dt=",1pE11.4)',it,time0,time

      end


************************************************************************
      subroutine DIFFUSION_IMPLICIT(time0,deltat,limiting,Nlim,Dlim,
     >                              verbose)
************************************************************************
      use GRID,ONLY: N=>Npoints,zz,d1l,d1m,d1r,d2l,d2m,d2r,
     >               dt_diff_im,xlower,xupper
      use PARAMETERS,ONLY: tfac,bc_low,bc_high,
     >                     outflux,inrate,outrate,vin,vout,
     >                     dust_diffuse
      use STRUCT,ONLY: Diff,nHtot,nHeps,crust_Neps,crust_Ncond
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use ELEMENTS,ONLY: NELEM,elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      implicit none
      real*8,intent(IN) :: time0,deltat
      logical,intent(in) :: limiting(NELEM)
      integer,intent(in) :: Nlim(NELEM),Dlim(NELEM,10),verbose
      integer :: i,j,it,e,el,Nstep,round
      real*8,dimension(N) :: xx,xnew,rest
      real*8,dimension(N,N) :: BB_1,BB_2,BB
      real*8 :: influx(NELEM)
      real*8 :: nD,time,dt
      character :: CR = CHAR(13)

      Nstep = 20
      time  = 0.d0
      dt    = deltat/Nstep
      call INIT_DIFFUSION(N,1,dt,BB_1)  ! for bc_low=1
      call INIT_DIFFUSION(N,2,dt,BB_2)  ! for bc_low=2
      
      do it=1,Nstep

        do round=1,2

          if (round==2) then
            !---------------------------------------------------------
            ! ***  figure out influx of the non-limiting elements  ***
            ! ***  from the influxes of the limiting elements and  ***
            ! ***  the stoichiometry of the crust                  ***
            !---------------------------------------------------------
            call INFLUXES(influx,limiting,Nlim,Dlim) 
          endif

          do e=1,NELM
            if (e==iel) cycle
            el = elnum(e)
            if (round==1.and..not.limiting(el)) cycle
            if (round==2.and.limiting(el)) cycle

            xx(:) = nHeps(el,:)/nHtot(:) 

            !------------------------------
            ! ***  boundary conditions  ***
            !------------------------------
            rest = xx
            nD = nHtot(1)*Diff(1)  
            if (limiting(el)) then   
              rest(1) = xlower(el)              ! fixed concentration
              bc_low = 1
              BB = BB_1
            else                              
              bc_low = 2
              rest(1) = -influx(el)/nD/d1l(1)   ! fixed flux
              BB = BB_2
            endif
            nD = nHtot(N)*Diff(N)  
            if (bc_high==1) then
              rest(N) = xupper(el) 
            else if (bc_high==2) then
              rest(N) = -outflux/nD/d1r(N)           
            else if (bc_high==3) then
              rest(N) = 0.d0
            endif

            !----------------------------
            ! ***  implicit timestep  ***
            !----------------------------
            do i=1,N
              xnew(i) = 0.d0 
              do j=1,N 
                xnew(i) = xnew(i) + BB(i,j)*rest(j)
              enddo  
            enddo  
            xx(:) = xnew(:)
            nHeps(el,:) = nHtot(:)*xnew(:)

            !-----------------------------------------------------------
            ! ***  compute fluxes through inner boundary in round=1  ***
            !-----------------------------------------------------------
            nD = nHtot(1)*Diff(1)  
            if (limiting(el)) then
              influx(el) = -nD
     >                * (d1l(1)*xx(1) + d1m(1)*xx(2) + d1r(1)*xx(3))
            endif
            !nD = nHtot(N)*Diff(N)  
            !if (bc_high==1) then
            !  outflux = -nD*(d1l(N)*xx(N-2) + d1m(N)*xx(N-1) + d1r(N)*xx(N)) 
            !else if (bc_high==2) then
            !else if (bc_high==3) then
            !  outflux = outrate*nHtot(N)*xx(N)*vout 
            !endif

          enddo  
        enddo  

        do e=1,NELM
          if (e==iel) cycle
          el = elnum(e)
          !----------------------------------------
          ! ***  update crust column densities  ***
          !----------------------------------------
          if (crust_Neps(el)>0.Q0) then   
            print*,elnam(el),REAL(crust_Neps(el)),influx(el)*dt
            crust_Neps(el) = crust_Neps(el) - influx(el)*dt
            if (crust_Neps(el)<0.Q0) then
              print*,elnam(el),"negative crust column density"
              stop
            endif  
          endif  
        enddo  

        time = time + dt
        if (verbose>0) write(*,'(TL10," DIFFUSION:",I8,A,$)') it,CR

      enddo
      if (verbose>0) print'(" DIFFUSION:",I8,"  time=",1pE11.4,
     >               "  Dt=",1pE11.4)',it,time0,time

      end


************************************************************************
      SUBROUTINE INFLUXES(influx,limiting,Nlim,Dlim) 
************************************************************************
      use STRUCT,ONLY: Diff,nHtot,nHeps,crust_Neps,crust_Ncond
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use ELEMENTS,ONLY: NELEM,elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      implicit none
      real*8,intent(inout) :: influx(NELEM)
      logical,intent(in) :: limiting(NELEM)
      integer,intent(in) :: Nlim(NELEM),Dlim(NELEM,10)
      integer :: i,j,k,e,el,el2,el3,s,Neq1,Neq2,ii,jj,NEXTRA
      integer :: extra_nel(NDUST),extra_Ns(NDUST),extra_s(NDUST)
      real*8 :: AA(NELM,NELM),sol(NELM),rhs(NELM)
      real*8 :: dNcond(NDUST),check(NELEM)     
      real*8 :: qual
      real*8 :: extra_el(NDUST,20),extra_nu(NDUST,20)
      character(len=500) :: extra_nam(NDUST)
      logical :: disconsider(NDUST),found

      disconsider(:) = .false.
      NEXTRA = 0
      do e=1,NELM
        if (e==iel) cycle
        el = elnum(e)
        if (Nlim(el).le.1) cycle
        NEXTRA = NEXTRA+1
        extra_nam(NEXTRA) = ''
        extra_nel(NEXTRA) = 0
        extra_Ns(NEXTRA) = 0
        do i=1,Nlim(el)
          s = Dlim(el,i)
          disconsider(s) = .true.
          extra_nam(NEXTRA) = trim(extra_nam(NEXTRA))
     >                        //"-"//trim(dust_nam(s))
          extra_Ns(NEXTRA) = extra_Ns(NEXTRA)+1
          extra_s(NEXTRA) = s
          do j=1,dust_nel(s)
            el2 = dust_el(s,j)
            found = .false.
            do k=1,extra_nel(NEXTRA)
              el3 = extra_el(NEXTRA,k)
              if (el2==el3) then
                found = .true.
                exit
              endif
            enddo
            if (found) then
              extra_nu(NEXTRA,k) = extra_nu(NEXTRA,k)+dust_nu(s,j)
            else
              extra_nel(NEXTRA) = extra_nel(NEXTRA)+1
              k = extra_nel(NEXTRA)
              extra_el(NEXTRA,k) = el2
              extra_nu(NEXTRA,k) = dust_nu(s,j)
            endif
          enddo
        enddo
      enddo
      do i=1,NEXTRA
        print*,trim(extra_nam(i))
        do j=1,extra_nel(i)
          el = extra_el(i,j)
          print*,elnam(el),int(extra_nu(i,j))
        enddo  
      enddo

      !--------------------------------------------
      ! ***  fill in equation system            ***
      ! ***  ((stoich))*(dNcond/dt) = (influx)  ***
      ! ***  for limiting element fluxes        ***
      !--------------------------------------------
      ii = 0
      do e=1,NELM
        if (e==iel) cycle
        el = elnum(e)
        if (.not.limiting(el)) cycle
        if (crust_Neps(el).le.0.Q0) cycle
        ii = ii+1
        Neq1 = ii
        rhs(ii) = influx(el)
        AA(ii,:) = 0.d0
        jj = 0
        do s=1,NDUST
          if (crust_Ncond(s)<=0.Q0) cycle
          if (disconsider(s)) cycle
          jj = jj+1
          if (ii==1) print*,jj,dust_nam(s)
          Neq2 = jj
          do i=1,dust_nel(s)
            el2 = dust_el(s,i)
            if (el==el2) AA(ii,jj) = dust_nu(s,i)
          enddo
        enddo
        do s=1,NEXTRA
          jj = jj+1
          if (ii==1) print*,jj,trim(extra_nam(s))
          Neq2 = jj
          do i=1,extra_nel(s)
            el2 = extra_el(s,i)
            if (el==el2) AA(ii,jj) = extra_nu(s,i)
          enddo
        enddo   
        print*,ii,elnam(el)
      enddo  
      print*,Neq1,Neq2
      do i=1,Neq1
        print'(99(I3))',int(AA(i,1:Neq2))
      enddo  
      if (Neq1.ne.Neq2) stop "Neq1<>Neq2 in diffusion."
      stop

      !--------------------------------
      ! ***  solve equation system  ***
      !--------------------------------
      call GAUSS8(NELM,Neq1,AA,sol,rhs)
      dNcond(:) = 0.d0
      jj = 0
      do s=1,NDUST
        if (crust_Ncond(s)<=0.Q0) cycle
        if (disconsider(s)) cycle
        jj = jj+1
        dNcond(s) = sol(jj)
      enddo  
      do i=1,NEXTRA
        jj = jj+1
        do j=1,extra_Ns(i)
          s = extra_s(i,j)
          dNcond(s) = sol(jj)
        enddo
      enddo  

      !----------------------------------------------------
      ! ***  determine fluxes of non-limiting elements  ***
      !----------------------------------------------------
      check = influx
      influx(:) = 0.d0
      do e=1,NELM
        if (e==iel) cycle
        el = elnum(e)
        if (crust_Neps(el).le.0.Q0) cycle
        do s=1,NDUST
          if (crust_Ncond(s)<=0.Q0) cycle
          do i=1,dust_nel(s)
            el2 = dust_el(s,i)
            if (el==el2) then
              influx(el) = influx(el)+dust_nu(s,i)*dNcond(s)
            endif  
          enddo
        enddo
        if (limiting(el)) then
          qual = ABS(influx(el)/check(el)-1.d0) 
          print'(A3,"  influx=",1pE12.4,"  check=",1pE9.2)',
     >                 elnam(el),influx(el),qual
        else
          print'(A3,"  influx=",1pE12.4)',elnam(el),influx(el)
        endif  
      enddo  

      end
