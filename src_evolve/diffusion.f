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
      real*8 :: nH,Tg,flux,fmin,dz,bsum,Ncol,rsort(NELEM)
      integer :: i,j,e,el,emin,ecount(NELEM),esort(NELEM)
      integer :: isort,ipass,Ncrust  
      logical :: in_crust(NELEM),limiting(NELEM)

      xlower = crust_gaseps
      eps(:) = nHeps(:,0)/nHtot(0)

      !------------------------------------------
      ! ***  identify most abundant elements  ***
      !------------------------------------------
      ecount(:) = 0
      do i=1,NDUST
        if (crust_Ncond(i)<=0.Q0) cycle
        do j=1,dust_nel(i)
          el = dust_el(i,j)
          ecount(el) = ecount(el)+1
        enddo  
      enddo
      
      esort(:) = 0
      rsort(:) = 0.d0
      isort = 0
      ipass = 0
      Ncrust = 0
      do e=1,NELM
        if (e==iel) cycle
        el = elnum(e)
        in_crust(el) = (crust_Neps(el)>0.Q0)
        isort = isort+1
        if (.not.in_crust(el)) then
          ipass = ipass+1 
          esort(ipass+1:isort) = esort(ipass:isort-1)
          rsort(ipass+1:isort) = rsort(ipass:isort-1)
          esort(ipass) = el
          rsort(ipass) = eps(el)
          cycle
        endif
        Ncrust = Ncrust+1
        i = ipass+1
        do i=ipass+1,isort-1
          if (eps(el)<rsort(i)) exit
          if (ecount(el)==1) exit
        enddo   
        esort(i+1:isort) = esort(i:isort-1)
        rsort(i+1:isort) = rsort(i:isort-1)
        esort(i) = el
        rsort(i) = eps(el)
        if (ecount(el)==1) rsort(i)=1.E-99
      enddo  
      do i=1,NDUST
        if (crust_Ncond(i)>0.Q0) Ncrust=Ncrust-1
      enddo
      limiting(:) = .true.
      limiting(isort+1-Ncrust:isort) = .false.
      do i=1,isort
        el = esort(i) 
        print'(A3,2(L2),2(1pE10.3))',elnam(el),in_crust(el),
     >                         limiting(i),rsort(i),eps(el) 
      enddo

      if (implicit.and.deltat>30*dt_diff_ex) then
        if (verbose>0) then
          print*
          print*,"entering IMPLICIT DIFFUSION ..."
          print*,"==============================="
        endif  
        call DIFFUSION_IMPLICIT(time,deltat,isort,ipass,esort,
     >                          in_crust,limiting,verbose)
      else
        if (verbose>1) then
          print*
          print*,"entering EXPLICIT DIFFUSION ..."
          print*,"==============================="
        endif  
        call DIFFUSION_EXPLICIT(time,deltat,isort,ipass,esort,
     >                          in_crust,limiting,verbose)
      endif
      end


************************************************************************
      subroutine DIFFUSION_EXPLICIT(time0,deltat,isort,ipass,esort,
     >                              in_crust,limiting,verbose)
************************************************************************
      use NATURE,ONLY: pi
      use GRID,ONLY: N=>Npoints,zz,zweight,d1l,d1m,d1r,d2l,d2m,d2r,
     >               dt_diff_ex,xlower,xupper
      use PARAMETERS,ONLY: bc_low,bc_high,
     >                     outflux,inrate,outrate,vin,vout
      use STRUCT,ONLY: Diff,nHtot,nHeps,
     >                 crust_Neps,crust_Ncond,crust_gaseps
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use ELEMENTS,ONLY: NELEM,elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      implicit none
      real*8,intent(in) :: time0
      real*8,intent(inout) :: deltat
      logical,intent(in) :: in_crust(NELEM),limiting(NELEM)
      integer,intent(in) :: isort,ipass,esort(NELEM),verbose
      integer :: i,it,e,el,round
      real*8,dimension(-2:N) :: xx,xold,rate
      real*8 :: influx(NELEM)
      real*8 :: D,nD,d1,d2,d1nD,time,dt,dz,dNcol,xl,xm,xr
      character :: CR = CHAR(13)
      character(len=1) :: char1
      logical :: IS_NAN,exhausted

      time = 0.d0

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
        
        round = 1
        influx(:) = 0.d0
        do e=1,isort

          el = esort(e)
          if ((round==1).and.(.not.limiting(e))) then
            round = 2 
            call INFLUXES(influx,isort,ipass,esort,limiting) 
          endif  

          !-------------------------------------------
          ! ***  d/dt(nH*x) = d/dz(nH*Diff*dx/dz)  ***
          !-------------------------------------------
          xx(:) = nHeps(el,:)/nHtot(:) 
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
            dNcol = dNcol + nHtot(i)*(xx(i)-xold(i))*zweight(i)
          enddo
          nD = nHtot(0)*Diff(0)  
          if (in_crust(el).and.limiting(e)) then   
            xl = 0.5*(xx(-1)+xold(-1)) 
            xm = 0.5*(xx( 0)+xold( 0)) 
            xr = 0.5*(xx(+1)+xold(+1)) 
            influx(el) = -nD*(d1l(0)*xl + d1m(0)*xm + d1r(0)*xr)
            print'(A3,"  influx=",2(1pE12.4))',
     >           elnam(el),influx(el),dNcol/dt
            influx(el) = dNcol/dt
          else   
            nD = nHtot(0)*Diff(0)
            xx(-1) = (-influx(el)/nD - d1m(0)*xx(0) - d1r(0)*xx(1))
     >               /d1l(0)                   ! constant flux 
            !print'(A3,"  influx=",99(1pE12.4))',elnam(el),-nD
     >      !     *(d1l(0)*xx(-1) + d1m(0)*xx(0) + d1r(0)*xx(1))
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

        exhausted = .false.
        do e=1,isort
          el = esort(e)
          !----------------------------------------
          ! ***  update crust column densities  ***
          !----------------------------------------
          if (in_crust(el)) then   
            print'(A4,2(1pE14.6))',elnam(el),
     >            crust_Neps(el),influx(el)*dt
            crust_Neps(el) = crust_Neps(el) - influx(el)*dt
            if (crust_Neps(el)<0.Q0) then
              exhausted = .true. 
              print*,elnam(el),"*** negative crust column density"
            endif  
          endif  

        enddo  

        time = time + dt
        if (verbose>0) write(*,'(TL10," DIFFUSION:",I8,A,$)') it,CR
        if (time.ge.deltat) exit
        if (exhausted) then
          deltat = time
          exit
        endif  

      enddo  
      if (verbose>0) print'(" DIFFUSION:",I8,"  time=",1pE11.4,
     >               "  Dt=",1pE11.4)',it,time0,time

      end


************************************************************************
      subroutine DIFFUSION_IMPLICIT(time0,deltat,isort,ipass,esort,
     >                              in_crust,limiting,verbose)
************************************************************************
      use GRID,ONLY: N=>Npoints,zz,zweight,d1l,d1m,d1r,d2l,d2m,d2r,
     >               dd1l,dd1m,dd1r,dt_diff_im,xlower,xupper
      use PARAMETERS,ONLY: bc_low,bc_high,
     >                     outflux,inrate,outrate,vin,vout
      use STRUCT,ONLY: Diff,nHtot,nHeps,crust_Neps,crust_Ncond
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use ELEMENTS,ONLY: NELEM,elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      implicit none
      real*8,intent(IN) :: time0
      real*8,intent(INOUT) :: deltat
      logical,intent(in) :: limiting(NELEM),in_crust(NELEM)
      integer,intent(in) :: isort,ipass,esort(NELEM),verbose
      integer :: i,j,it,e,el,Nstep,round
      real*8,dimension(0:N) :: xx,xnew,rest
      real*8,dimension(N+1,N+1) :: BB_1,BB_2,BB
      real*8 :: nHold(NELEM,-2:N),influx(NELEM),dNcol,Natmos
      real*8 :: nD,time,dt,xl,xm,xr
      character :: CR = CHAR(13)
      character(len=1) :: char1
      logical :: exhausted,toomuch,reduced

      Nstep = 20
      time  = 0.d0
      dt    = deltat/Nstep
      reduced = .false.
      call INIT_DIFFUSION(N,1,dt,BB_1)  ! for bc_low=1
      call INIT_DIFFUSION(N,2,dt,BB_2)  ! for bc_low=2
      
      do it=1,Nstep

 100    continue
        round = 1
        influx(:) = 0.d0
        nHold(:,:) = nHeps(:,:)
        if (MINVAL(nHeps)<0.d0) then
          print*,"*** negative nHeps it=",it
          stop
        endif  

        do e=1,isort

          el = esort(e)
          if ((round==1).and.(.not.limiting(e))) then
            round = 2 
            call INFLUXES(influx,isort,ipass,esort,limiting) 
          endif  
          xx(0:N) = nHeps(el,0:N)/nHtot(0:N) 

          !------------------------------
          ! ***  boundary conditions  ***
          !------------------------------
          rest = xx
          nD = nHtot(0)*Diff(0)  
          if (in_crust(el).and.limiting(e)) then   
            rest(0) = xlower(el)                 ! fixed concentration
            bc_low = 1
            BB = BB_1
          else                              
            bc_low = 2
            rest(0) = -influx(el)/nD/dd1l(0)     ! fixed flux
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
          do i=0,N
            xnew(i) = 0.d0 
            do j=0,N 
              xnew(i) = xnew(i) + BB(i+1,j+1)*rest(j)
            enddo  
          enddo  
          print'(A3,99(1pE10.3))',elnam(el),xx(0:3)
          print'(3x,99(1pE10.3))',xnew(0:3)

          !------------------------------------------------
          ! ***  compute fluxes through inner boundary  ***
          !------------------------------------------------
          dNcol = 0.d0
          do i=0,N
            dNcol = dNcol + nHtot(i)*(xnew(i)-xx(i))*zweight(i)
          enddo

          nD = nHtot(0)*Diff(0)  
          if (in_crust(el).and.limiting(e)) then   
            xl = 0.5*(xx(0)+xnew(0)) 
            xm = 0.5*(xx(1)+xnew(1)) 
            xr = 0.5*(xx(2)+xnew(2)) 
            influx(el) = -nD*(dd1l(0)*xl + dd1m(0)*xm + dd1r(0)*xr)
            print'(A3,"  influx=",2(1pE12.4))',
     >           elnam(el),influx(el),dNcol/dt
            influx(el) = dNcol/dt
          endif  

          xx(:) = xnew(:)
          nHeps(el,0:N) = nHtot(0:N)*xnew(0:N)

        enddo  
        
        !----------------------------------------------------
        ! ***  check whether influx would substantially   *** 
        ! ***  increase total atmospheric column density  ***
        !----------------------------------------------------
        toomuch = .false.
        do e=1,isort
          el = esort(e)
          if (in_crust(el)) then   
            Natmos = 0.d0
            do i=0,N
              Natmos = Natmos + nHold(el,i)*zweight(i)
            enddo
            if (influx(el)*dt>0.2*Natmos) then
              print*,"*** "//elnam(el)//" too large timestep",
     >               influx(el)*dt/Natmos
              dt = dt/2.0
              toomuch = .true.
              exit
            endif  
          endif
        enddo     
        if (toomuch) then
          nHeps(:,:) = nHold(:,:) 
          call INIT_DIFFUSION(N,1,dt,BB_1)
          call INIT_DIFFUSION(N,2,dt,BB_2)
          read'(A1)',char1
          reduced = .true.
          goto 100
        endif  

        exhausted = .false.
        do e=1,NELM
          if (e==iel) cycle
          el = elnum(e)
          !----------------------------------------
          ! ***  update crust column densities  ***
          !----------------------------------------
          if (in_crust(el)) then 
            print*,elnam(el),REAL(crust_Neps(el)),influx(el)*dt
            crust_Neps(el) = crust_Neps(el) - influx(el)*dt
            if (crust_Neps(el)<0.Q0) then
              exhausted = .true. 
              print*,elnam(el),"*** negative crust column density"
              read'(A1)',char1
            endif  
          endif  
        enddo  

        time = time + dt
        if (verbose>0) write(*,'(TL10," DIFFUSION:",I8,A,$)') it,CR
        if (exhausted) exit

      enddo
      if (verbose>0) print'(" DIFFUSION:",I8,"  time=",1pE11.4,
     >               "  Dt=",1pE11.4)',it,time0,time

      deltat = time    ! actually advanced timestep

      end


************************************************************************
      SUBROUTINE INFLUXES(influx,isort,ipass,esort,limiting) 
************************************************************************
      use STRUCT,ONLY: crust_Neps,crust_Ncond
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use ELEMENTS,ONLY: NELEM,elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      implicit none
      real*8,intent(inout) :: influx(NELEM)
      logical,intent(in) :: limiting(NELEM)
      integer,intent(in) :: isort,ipass,esort(NELEM)
      integer :: i,ii,jj,e,el,el2,i1,i2,s,Neq1,Neq2
      real*8 :: AA(NELM,NELM),sol(NELM),rhs(NELM)
      real*8 :: dNcond(NDUST),check(NDUST),qual

      !--------------------------------------------
      ! ***  fill in equation system            ***
      ! ***  ((stoich))*(dNcond/dt) = (influx)  ***
      ! ***  for limiting element fluxes        ***
      !--------------------------------------------
      i1 = ipass+1
      do i=i1,isort
        if (limiting(i)) i2=i 
      enddo 
      Neq1 = i2-i1+1

      ii = 0
      do e=i1,i2
        el = esort(e)
        ii = ii+1
        rhs(ii) = influx(el)
        AA(ii,:) = 0.d0
        jj = 0
        do s=1,NDUST
          if (crust_Ncond(s)<=0.Q0) cycle
          jj = jj+1
          if (e==i1) print*,jj,dust_nam(s)
          Neq2 = jj
          do i=1,dust_nel(s)
            el2 = dust_el(s,i)
            if (el==el2) AA(ii,jj) = dust_nu(s,i)
          enddo
        enddo
        print*,ii,elnam(el)
      enddo  
      print*,Neq1,Neq2
      do i=1,Neq1
        print'(99(I3))',int(AA(i,1:Neq2))
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
        jj = jj+1
        dNcond(s) = sol(jj)
      enddo  

      !----------------------------------------------------
      ! ***  determine fluxes of non-limiting elements  ***
      !----------------------------------------------------
      check = influx
      influx(:) = 0.d0
      do e=i1,isort
        el = esort(e)
        do s=1,NDUST
          if (crust_Ncond(s)<=0.Q0) cycle
          do i=1,dust_nel(s)
            el2 = dust_el(s,i)
            if (el==el2) then
              influx(el) = influx(el)+dust_nu(s,i)*dNcond(s)
            endif  
          enddo
        enddo
        if (limiting(e)) then
          qual = ABS(influx(el)/check(el)-1.d0) 
          print'(A3,"  influx=",1pE12.4,"  check=",1pE9.2)',
     >                 elnam(el),influx(el),qual
        else
          print'(A3,"  influx=",1pE12.4)',elnam(el),influx(el)
        endif  
      enddo  

      end
