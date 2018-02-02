************************************************************************
      subroutine DIFFUSION(time,deltat,verbose)
************************************************************************
      use PARAMETERS,ONLY: implicit
      use GRID,ONLY: zz,xlower,xupper,dt_diff_ex
      use STRUCT,ONLY: Temp,Diff,nHtot,nHeps,
     >                 crust_Neps,crust_Ncond,crust_gaseps
      use ELEMENTS,ONLY: NELEM,elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use EXCHANGE,ONLY: nat,nmol
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in) :: time
      real*8,intent(inout) :: deltat
      integer,intent(in) :: verbose
      real(kind=qp) :: eps(NELEM),Sat(NDUST)
      real*8 :: nH,Tg,flux,fmin,dz,bsum,Ncol,fac=0.05
      real*8 :: branch(NELEM,10)
      integer :: s,i,j,e,el,el2,emin
      integer :: Nlim(NELEM),Dlim(NELEM,10)
      logical :: evaporates(NDUST),limiting(NELEM)

      xlower = crust_gaseps
      eps(:) = nHeps(:,0)/nHtot(0)
      dz = zz(1)-zz(0)
      do e=1,NELM
        if (e==iel) cycle
        el = elnum(e)
        !if (crust_Neps(el)>0.d0) cycle
        flux = -Diff(0)*nHtot(0)*(eps(el)-crust_gaseps(el))/dz
        Ncol = nHeps(el,0)*dz + crust_Neps(el)
        if (ABS(flux)*deltat > fac*Ncol) then
          deltat = fac*Ncol/ABS(flux)
          print'(A3," coldens=",1pE10.3," change=",1pE10.3,
     >           "  dt=",1pE10.3)',elnam(el),Ncol,
     >           ABS(flux)*deltat/Ncol,deltat
        else
          print'(A3," coldens=",1pE10.3," change=",1pE10.3)',
     >           elnam(el),Ncol,ABS(flux)*deltat/Ncol
        endif
      enddo  

      if (.false.) then
      !------------------------------------------------
      ! ***  Consider layer=2 to determine whether  ***
      ! ***  solids grow or evaporate.              ***
      ! ***  Determine flux limiting elements.      ***
      !------------------------------------------------
      Tg = Temp(1)
      nH = nHtot(1)
      call GGCHEM(nH,Tg,eps,.false.,0)
      call SUPERSAT(Tg,nat,nmol,Sat)
      evaporates(:) = .false.
      limiting(:) = .false.
      Nlim(:) = 0
      do s=1,NDUST
        if (crust_Ncond(s)<=0.Q0) cycle
        if (Sat(s)<1.Q0) evaporates(s)=.true. 
        print'(A15,"  evaporates=",L1)',dust_nam(s),evaporates(s)
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
      !-------------------------------------------------------
      ! ***  if one element limites the growth of several  ***
      ! ***  materials, fix branching ratios dependent on  ***
      ! ***  the abundance of the second most important    ***
      ! ***  constituent element                           ***
      !-------------------------------------------------------
      do el=1,NELEM
        if (Nlim(el).le.1) cycle
        bsum = 0.d0
        do i=1,Nlim(el)
          branch(el,i) = 1.d0
          s = Dlim(el,i)
          do j=1,dust_nel(s)
            el2 = dust_el(s,j)
            if (el2==el) cycle
            branch(el,i) = MIN(branch(el,i),
     >                         crust_gaseps(el2)/dust_nu(s,j))
          enddo
          bsum = bsum + branch(el,i)
        enddo
        branch(el,:) = branch(el,:)/bsum
        print*,elnam(el),": ",dust_nam(Dlim(el,1:Nlim(el)))
        print*,branch(el,1:Nlim(el))
      enddo 
      endif

      if (implicit.and.deltat>100*dt_diff_ex) then
        if (verbose>1) then
          print*
          print*,"entering IMPLICIT DIFFUSION ..."
          print*,"==============================="
        endif  
        stop
        call DIFFUSION_IMPLICIT(time,deltat,limiting,Nlim,Dlim,
     >                          branch,verbose)
      else
        if (verbose>1) then
          print*
          print*,"entering EXPLICIT DIFFUSION ..."
          print*,"==============================="
        endif  
        call DIFFUSION_EXPLICIT(time,deltat,limiting,Nlim,Dlim,
     >                          branch,verbose)
      endif
      end


************************************************************************
      subroutine DIFFUSION_EXPLICIT(time0,deltat,limiting,Nlim,Dlim,
     >                              branch,verbose)
************************************************************************
      use NATURE,ONLY: pi
      use GRID,ONLY: N=>Npoints,zz,d1l,d1m,d1r,d2l,d2m,d2r,dt_diff_ex,
     >               xlower,xupper
      use PARAMETERS,ONLY: dust_diffuse,bc_low,bc_high,
     >                     outflux,inrate,outrate,vin,vout
      use STRUCT,ONLY: Diff,nHtot,nHeps,
     >                 crust_Neps,crust_Ncond,crust_gaseps
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use ELEMENTS,ONLY: NELEM,elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      implicit none
      real*8,intent(in) :: time0,deltat,branch(NELEM,10)
      logical,intent(in) :: limiting(NELEM)
      integer,intent(in) :: Nlim(NELEM),Dlim(NELEM,10),verbose
      integer :: i,it,e,el,round
      real*8,dimension(-2:N) :: xx,xold,rate
      real*8 :: influx(NELEM)
      real*8 :: D,nD,d1,d2,d1nD,time,dt,dz,dNcol
      character :: CR = CHAR(13)
      character(len=1) :: char1
      logical :: in_crust(NELEM),IS_NAN

      time = 0.d0
      do e=1,NELM
        if (e==iel) cycle
        el = elnum(e)
        in_crust(el) = (crust_Neps(el)>0.Q0)
      enddo  

      do it=1,9999999

        dt = MIN(dt_diff_ex,deltat-time)

        !------------------------------
        ! ***  boundary conditions  ***  
        !------------------------------
        do e=1,NELM
          if (e==iel) cycle
          el = elnum(e)
          if (in_crust(el)) then
            !--- put crust_eps on guard cells ---
            nHeps(el,-2) = nHtot(-2)*crust_gaseps(el)
            nHeps(el,-1) = nHtot(-1)*crust_gaseps(el)
          else
            !--- mirror structure on guard cells to enforce flux=0 ---
            nHeps(el,-2) = nHeps(el,2)/nHtot(2)*nHtot(-2)
            nHeps(el,-1) = nHeps(el,1)/nHtot(1)*nHtot(-1)
          endif
        enddo  

        
        do e=1,NELM
          if (e==iel) cycle
          el = elnum(e)
          !if (round==1.and..not.limiting(el)) cycle
          !if (round==2.and.limiting(el)) cycle

          xx(:) = nHeps(el,:)/nHtot(:) 

          !-------------------------------------------
          ! ***  d/dt(nH*x) = d/dz(nH*Diff*dx/dz)  ***
          !-------------------------------------------
          rate(:) = 0.0
          do i=-1,N-1
            nD   = nHtot(i)*Diff(i)  
            d1   = d1l(i)*xx(i-1) + d1m(i)*xx(i) + d1r(i)*xx(i+1)
            d2   = d2l(i)*xx(i-1) + d2m(i)*xx(i) + d2r(i)*xx(i+1)
            d1nD = d1l(i)*nHtot(i-1)*Diff(i-1) 
     >           + d1m(i)*nHtot(i)  *Diff(i) 
     >           + d1r(i)*nHtot(i+1)*Diff(i+1) 
            rate(i) = nD*d2 + d1nD*d1
          enddo
          
          !----------------------------
          ! ***  explicit timestep  ***
          !----------------------------
          xold = xx
          do i=-1,N-1
            xx(i) = xx(i) + rate(i)*dt/nHtot(i)
          enddo  

          !------------------------------
          ! ***  boundary conditions  ***
          !------------------------------
          dNcol = 0.d0
          do i=0,N
            dz = zz(i)-zz(i-1) 
            dNcol = dNcol + nHtot(i)*(xx(i)-xold(i))*dz
          enddo
          nD = nHtot(0)*Diff(0)  
          if (in_crust(el)) then   
            influx(el) = -nD
     >              *(d1l(0)*xx(-1) + d1m(0)*xx(0) + d1r(0)*xx(1))
            print'(A3,"  influx=",2(1pE12.4))',
     >           elnam(el),influx(el),dNcol/dt
            influx(el) = dNcol/dt
          else   
            influx(el) = 0.d0 
            nD = nHtot(0)*Diff(0)
            xx(-1) = (-influx(el)/nD - d1m(0)*xx(0) - d1r(0)*xx(1))
     >               /d1l(0)                   ! constant flux 
            print'(A3,"  influx=",99(1pE12.4))',elnam(el),-nD
     >           *(d1l(0)*xx(-1) + d1m(0)*xx(0) + d1r(0)*xx(1))
          endif  
          nD = nHtot(N)*Diff(N)  
          if (bc_high==1) then
            xx(N) = xupper(el)                ! const concentration
            outflux = -nD
     >           *(d1l(N)*xx(N-2) + d1m(N)*xx(N-1) + d1r(N)*xx(N)) 
          else if (bc_high==2) then   
            xx(N) = (-outflux/nD - d1l(N)*xx(N-2) - d1m(N)*xx(N-1))
     >            /d1r(N)                     ! constant flux 
          else if (bc_high==3) then   
            outflux = nHtot(N)*xx(N)*outrate*vout  
            xx(N) = -(d1l(N)*xx(N-2) + d1m(N)*xx(N-1))
     >           /(d1r(N) + outrate*vout/Diff(N))
          endif   

          !------------------------------------------
          ! ***  map solution on atmosphere grid  ***
          !------------------------------------------
          print'(A3,7(1pE10.3))',elnam(el),xold(-2:3)
          print'(3x,7(1pE10.3))',xx(-2:3)
          do i=-2,N
            nHeps(el,i) = nHtot(i)*xx(i)
            if (xx(i)<0.d0) then
              print*,"negative elem.abund. in diffusion." 
              print*,elnam(el),i
              read'(A1)',char1
              if (i==N) then 
                nHeps(el,i) = nHtot(i)*xx(N-1)
              else  
                nHeps(el,i) = nHtot(i)*1.E-50
              endif
            endif    
          enddo 

        enddo  

        do e=1,NELM
          if (e==iel) cycle
          el = elnum(e)
          !----------------------------------------
          ! ***  update crust column densities  ***
          !----------------------------------------
          if (in_crust(el)) then   
            print'(A4,2(1pE14.6))',elnam(el),
     >            crust_Neps(el),influx(el)*dt
            crust_Neps(el) = crust_Neps(el) - influx(el)*dt
            if (crust_Neps(el)<0.Q0) then
              print*,elnam(el),"*** negative crust column density"
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
     >                              branch,verbose)
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
      real*8,intent(IN) :: time0,deltat,branch(NELEM,10)
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
            call INFLUXES(influx,limiting,Nlim,Dlim,branch) 
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
            print'(A3,4(1pE10.3))',elnam(el),xx(1:4)
            print'(3x,4(1pE10.3))',xnew(1:4)
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
      SUBROUTINE INFLUXES(influx,limiting,Nlim,Dlim,branch) 
************************************************************************
      use STRUCT,ONLY: Diff,nHtot,nHeps,crust_Neps,crust_Ncond
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use ELEMENTS,ONLY: NELEM,elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      implicit none
      real*8,intent(inout) :: influx(NELEM)
      logical,intent(in) :: limiting(NELEM)
      integer,intent(in) :: Nlim(NELEM),Dlim(NELEM,10)
      real*8,intent(in) :: branch(NELEM,10)
      integer :: i,j,k,e,el,el2,el3,s,Neq1,Neq2,ii,jj,NEXTRA
      integer :: extra_nel(NDUST),extra_Ns(NDUST),extra_s(NDUST,20)
      real*8 :: AA(NELM,NELM),sol(NELM),rhs(NELM)
      real*8 :: dNcond(NDUST),check(NELEM)     
      real*8 :: qual
      real*8 :: extra_el(NDUST,20),extra_nu(NDUST,20)
      character(len=500) :: extra_nam(NDUST)
      logical :: disconsider(NDUST),found

      do e=1,NELM
        if (e==iel) cycle
        el = elnum(e)
        if (limiting(el)) cycle
        influx(el) = 0.d0
      enddo
      return

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
          extra_s(NEXTRA,i) = s
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
              extra_nu(NEXTRA,k) = extra_nu(NEXTRA,k)
     >                           + branch(el,i)*dust_nu(s,j)
            else
              extra_nel(NEXTRA)  = extra_nel(NEXTRA)+1
              k = extra_nel(NEXTRA)
              extra_el(NEXTRA,k) = el2
              extra_nu(NEXTRA,k) = branch(el,i)*dust_nu(s,j)
            endif
          enddo
        enddo
      enddo
      do i=1,NEXTRA
        print*,trim(extra_nam(i))
        do j=1,extra_nel(i)
          el = extra_el(i,j)
          print*,elnam(el),extra_nu(i,j)
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
        print'(99(1pE10.3))',AA(i,1:Neq2)
      enddo  
      if (Neq1.ne.Neq2) stop "Neq1<>Neq2 in diffusion."

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
      do e=1,NELM
        if (e==iel) cycle
        el = elnum(e)
        if (Nlim(el).le.1) cycle
        jj = jj+1
        do i=1,Nlim(el)
          s = Dlim(el,i)
          dNcond(s) = branch(el,i)*sol(jj)
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
