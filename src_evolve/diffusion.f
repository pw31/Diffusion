************************************************************************
      subroutine DIFFUSION(time,deltat,reduced,Nsolve,indep)
************************************************************************
      use PARAMETERS,ONLY: implicit,Rplanet,logg,verbose,bc_high
      use GRID,ONLY: zz,xlower,xupper,dt_diff_ex,Npoints
      use STRUCT,ONLY: Temp,Diff,nHtot,nHeps,
     >                 crust_Neps,crust_Ncond,crust_gaseps
      use ELEMENTS,ONLY: NELEM,elnam,eps0
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el,molmass,NMOLE,cmol
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use EXCHANGE,ONLY: nat,nmol,nel,Fe,O,H
      use NATURE,ONLY: bk,pi
      use JEANS_ESCAPE,ONLY: EXTRA,Hfrac,jpern
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in) :: time
      real*8,intent(inout) :: deltat
      logical,intent(out) :: reduced
      integer,intent(in) :: Nsolve,indep(NELEM)
      real(kind=qp) :: eps(NELEM),Sat(NDUST)
      real*8 :: nH,Tg,flux,fmin,dz,bsum,Ncol,rsort(NELEM+EXTRA),sort
      real*8 :: rtmp,abun
      integer :: i,j,e,el,emin,esort(NELEM+EXTRA)
      integer :: epure(NELEM),ecount(NELEM),epos(NELEM)
      integer :: isort,ipass,Ncond,Ndep,i1,i2,e1,e2,itmp  
      logical :: reduced1,reduced2
      logical :: in_crust(NELEM+EXTRA),limiting(NELEM+EXTRA)
      logical :: act(NDUST),elim(NELEM),eflag(NELEM)
      logical,save :: firstCall=.true.
      integer,save :: iFe_l=0,iFeO_l=0
      character(len=3),save :: spnam(NELEM+EXTRA)

      if (firstCall) then
        do i=1,NDUST
          if (dust_nam(i).eq.'Fe[l]')  iFe_l=i 
          if (dust_nam(i).eq.'FeO[l]') iFeO_l=i 
        enddo
        spnam(1:NELEM) = elnam(1:NELEM)
        spnam(NELEM+1) = 'Hat'
        spnam(NELEM+2) = 'H2'
        spnam(NELEM+3) = 'H2O'
        firstCall = .false.
      endif
  
      xlower = crust_gaseps
      eps(:) = nHeps(:,0)/nHtot(0)

      !-------------------------------------------------------------
      ! ***  identify most abundant elements for inner boundary  ***
      !-------------------------------------------------------------
      epure(:)  = 0
      ecount(:) = 0
      act(:) = .false.
      do i=1,NDUST
        if (crust_Ncond(i)<=0.Q0) cycle
        act(i) = .true.
        if (dust_nel(i)==1) epure(dust_el(i,1))=1
        do j=1,dust_nel(i)
          el = dust_el(i,j)
          ecount(el) = ecount(el)+1
        enddo  
      enddo
      esort(:) = 0
      rsort(:) = 0.d0
      isort = 0
      ipass = 0
      Ndep  = 0
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
        Ndep = Ndep+1
        i = ipass+1
        sort = eps(el)
        if (epure(el)==1) sort=1.E-99   ! keep the pure atomic species
        do i=ipass+1,isort-1
          if (sort<rsort(i)) exit
        enddo   
        esort(i+1:isort) = esort(i:isort-1)
        rsort(i+1:isort) = rsort(i:isort-1)
        esort(i) = el
        rsort(i) = sort
      enddo
      Ncond = 0
      do i=1,NDUST
        if (crust_Ncond(i)>0.Q0) then
          Ncond = Ncond+1
          Ndep = Ndep-1
        endif  
      enddo
      limiting(:) = .true.
      limiting(isort+1-Ndep:isort) = .false.
      elim(:) = .true.
      do i=1,isort
        el = esort(i) 
        epos(el) = i
        if (i>=isort+1-Ndep) elim(el)=.false.
      enddo
      !-----  correction special cases ...  ------
      if (act(iFe_l).and.act(iFeO_l).and.ecount(Fe)==2.and.
     >    (.not.elim(O))) then
        !--- O cannot be non-limiting element---
        i1 = epos(O)
        i2 = isort-Ndep 
        e1 = esort(i1)
        e2 = esort(i2)
        print*,"diffusion: swapping "//elnam(e1)//" with "//elnam(e2)
        itmp = esort(i1)
        rtmp = rsort(i1)
        esort(i1) = esort(i2)
        rsort(i1) = rsort(i2)
        esort(i2) = itmp
        rsort(i2) = rtmp
        elim(e1)  = .true.
        elim(e2)  = .false.        
      endif   
      !-----  better use sorting from equil_cond if possible -----
      print*,"sorting and type of elements ..."
      if (Ncond==Nsolve) then
        print*,"before: ",elnam(esort(1:isort)) 
        !print*,elnam(indep(1:Nsolve)) 
        !print*,isort-Ncond-Ndep,Ncond,Ndep
        eflag = .false.
        eflag(esort(1:isort-Ncond-Ndep)) = .true.
        eflag(indep(1:Nsolve)) = .true.
        esort(isort-Ncond-Ndep+1:isort-Ndep) = indep(1:Nsolve)
        j=isort-Ndep
        do e=1,NELM
          if (e==iel) cycle
          el = elnum(e)
          if (eflag(el)) cycle
          j=j+1
          esort(j) = el
        enddo  
        print*,"after:  ",elnam(esort(1:isort)) 
      endif
   
      !--- extra slots for detailed H ---
      do e=1,isort
        el = esort(e)
        if (el==H) exit
      enddo
      esort(e+1+EXTRA:isort+EXTRA) = esort(e+1:isort)
      esort(e+1) = NELEM+1
      esort(e+2) = NELEM+2
      esort(e+3) = NELEM+3
      in_crust(NELEM+1:NELEM+3) = in_crust(H)
      limiting(e+1+EXTRA:isort+EXTRA) = limiting(e+1:isort)
      limiting(e+1) = limiting(e)
      limiting(e+2) = limiting(e)
      limiting(e+3) = limiting(e)
      if (.not.in_crust(H)) ipass=ipass+EXTRA
      isort = isort+EXTRA

      !---------------------------------------
      ! ***  upper boundary: Jeans Escape  ***
      !---------------------------------------
      if (bc_high==3) then
        call ESCAPE(deltat,reduced1,0)
      endif  

      if (verbose>=0) then
        do i=1,isort
          el = esort(i) 
          if (el<=NELEM)   abun=eps(el)
          if (el==NELEM+1) abun=eps(H)*Hfrac(1)
          if (el==NELEM+2) abun=eps(H)*Hfrac(2)
          if (el==NELEM+3) abun=eps(H)*Hfrac(3)
          print'(A4,2(L3),99(1pE11.3))',spnam(el),in_crust(el)
     >                             ,limiting(i),abun,jpern(el)
        enddo
      endif  

      if (implicit.and.deltat>30*dt_diff_ex) then
        if (verbose>=0) then
          print*
          print*,"entering IMPLICIT DIFFUSION ..."
          print*,"==============================="
        endif  
        call DIFFUSION_IMPLICIT(time,deltat,isort,ipass,esort,
     >                          in_crust,limiting,reduced2,spnam)
      else
        if (verbose>=0) then
          print*
          print*,"entering EXPLICIT DIFFUSION ..."
          print*,"==============================="
        endif  
        call DIFFUSION_EXPLICIT(time,deltat,isort,ipass,esort,
     >                          in_crust,limiting,reduced2,spnam)
      endif
      reduced = (reduced1.or.reduced2)

      end


************************************************************************
      subroutine DIFFUSION_EXPLICIT(time0,deltat,isort,ipass,esort,
     >                              in_crust,limiting,reduced,spnam)
************************************************************************
      use NATURE,ONLY: pi
      use GRID,ONLY: N=>Npoints,zz,zweight,d1l,d1m,d1r,d2l,d2m,d2r,
     >               dt_diff_ex,xlower,xupper
      use PARAMETERS,ONLY: bc_low,bc_high,verbose,dt_increase,
     >                     outflux,detailed_H
      use STRUCT,ONLY: Diff,nHtot,nHeps,
     >                 crust_Neps,crust_Ncond,crust_gaseps
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use ELEMENTS,ONLY: NELEM,elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      use JEANS_ESCAPE, ONLY: EXTRA,jpern,Hfrac
      use EXCHANGE,ONLY: H
      implicit none
      real*8,intent(in) :: time0
      real*8,intent(inout) :: deltat
      logical,intent(in)  :: in_crust(NELEM+EXTRA),limiting(NELEM+EXTRA)
      integer,intent(in)  :: isort,ipass,esort(NELEM+EXTRA)
      logical,intent(out) :: reduced
      character(len=3),intent(in) :: spnam(NELEM+EXTRA)
      integer :: i,it,e,el,round
      real*8,dimension(-2:N) :: xx,xold,xHtot,nHmerk,rate
      real*8 :: influx(NELEM+EXTRA),Hinflux,Natmos
      real*8 :: D,nD,d1,d2,d1nD,time,dt,dz,dNcol,xl,xm,xr
      real*8 :: nHold(NELEM+EXTRA,-2:N),crust_Nold(NELEM)
      character :: CR = CHAR(13)
      character(len=1) :: char1
      logical :: IS_NAN,exhausted,toomuch

      reduced = .false.
      exhausted = .false.
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
        nHold(:,:) = nHeps(:,:)
        crust_Nold(:) = crust_Neps(:)

 100    continue
        print*,"dt=",dt
        round = 1
        xHtot(:) = 0.d0
        influx(:) = 0.d0
        Hinflux = 0.d0
        toomuch = .false.

        do e=1,isort

          el = esort(e)
          if ((round==1).and.(.not.limiting(e))) then
            round = 2 
            call INFLUXES(influx,isort,ipass,esort,limiting) 
          endif  

          !-------------------------------------------
          ! ***  d/dt(nH*x) = d/dz(nH*Diff*dx/dz)  ***
          !-------------------------------------------
          if (el<=NELEM) then
            xx(:) = nHeps(el,:)/nHtot(:)                ! regular element
          else  
            xx(:) = nHeps(H,:)/nHtot(:)*Hfrac(el-NELEM) ! H,H2,H2O
          endif  
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

          if (round==2) then
            print*,"deviation "//spnam(el),xold(0)/xx(0)-1.d0
            if (ABS(xold(0)/xx(0)-1.d0)>0.01) then
              print*,'*** too large deviation'
              toomuch = .true.
            endif  
          endif

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
     >           spnam(el),influx(el),dNcol/dt
            influx(el) = dNcol/dt
          else   
            nD = nHtot(0)*Diff(0)
            xx(-1) = (-influx(el)/nD - d1m(0)*xx(0) - d1r(0)*xx(1))
     >               /d1l(0)                   ! constant flux 
            if (in_crust(el)) print'(A3,"  influx=",2(1pE12.4))',
     >           spnam(el),influx(el)
          endif  
          nD = nHtot(N)*Diff(N)  
          if (bc_high==1) then
            xx(N) = xupper(el)                    ! const concentration
            outflux = -nD
     >           *(d1l(N)*xx(N-2) + d1m(N)*xx(N-1) + d1r(N)*xx(N)) 
          else if (bc_high==2) then   
            xx(N) = (-outflux/nD - d1l(N)*xx(N-2) - d1m(N)*xx(N-1))
     >            /d1r(N)                         ! constant flux 
          else if (bc_high==3) then 
            outflux = nHtot(N)*xx(N)*jpern(el)    ! constant escape rate
            xx(N) = -(d1l(N)*xx(N-2) + d1m(N)*xx(N-1))
     >           /(d1r(N) + jpern(el)/Diff(N))
          endif   

          !------------------------------------------
          ! ***  map solution on atmosphere grid  ***
          !------------------------------------------
          if (verbose>0) then
            print'(A3,7(1pE10.3))',spnam(el),xold(-2:3)
            print'(3x,7(1pE10.3))',xx(-2:3)
          endif  
          do i=-2,N
            if (xx(i)<0.d0) then
              print*,"*** negative elem.abund. in diffusion." 
              print*,spnam(el),i
              if (verbose>0) read'(A1)',char1
              if (i==N) then 
                xx(i) = xx(N-1)
              else  
                xx(i) = 1.E-50
              endif
            endif    
          enddo
          if (el==H) then
            nHmerk(:) = nHtot(:)*xx(:)          ! the result for mean hydrogen
            print'("outflux <H>=",3(1pE15.7))',
     >           nHmerk(N),jpern(H),nHmerk(N)*jpern(H)
          else if (el<=NELEM) then
            nHeps(el,:) = nHtot(:)*xx(:)
          else
            nHeps(el,:) = nHtot(:)*xx(:)
            xHtot(:) = xHtot(:) + xx(:)         ! add together detailed hydrogen
            Hinflux  = Hinflux  + influx(el)
          endif  
          if (el==NELEM+EXTRA) then
            if (detailed_H) then 
              !print*,nHtot(:)*xHtot(:)/nHmerk(:) 
              nHeps(H,:) = nHtot(:)*xHtot(:)
              influx(H)  = Hinflux
              outflux = SUM(nHeps(NELEM+1:NELEM+EXTRA,N)
     >                     *jpern(NELEM+1:NELEM+EXTRA))
              print'("outflux Hat=",3(1pE15.7))',nHeps(NELEM+1,N),
     >             jpern(NELEM+1),nHeps(NELEM+1,N)*jpern(NELEM+1)
              print'("outflux H2 =",3(1pE15.7))',nHeps(NELEM+2,N),
     >             jpern(NELEM+2),nHeps(NELEM+2,N)*jpern(NELEM+2)
              print'("outflux H2O=",3(1pE15.7))',nHeps(NELEM+3,N),
     >             jpern(NELEM+3),nHeps(NELEM+3,N)*jpern(NELEM+3)
              print'("outflux <H>=",3(1pE15.7))',nHeps(H,N),outflux
            else  
              nHeps(H,:) = nHmerk(:)
            endif  
          endif  

        enddo  

        do e=1,isort
          el = esort(e)
          if (el>NELEM) cycle
          !----------------------------------------
          ! ***  update crust column densities  ***
          !----------------------------------------
          if (in_crust(el)) then   
            if (verbose>0) then 
              print'(A4,2(1pE14.6))',elnam(el),
     >              crust_Neps(el),influx(el)*dt
            endif  
            crust_Neps(el) = crust_Neps(el) - influx(el)*dt
            if (crust_Neps(el)<0.Q0) then
              exhausted = .true. 
              reduced = .true.
              print*,elnam(el),"*** negative crust column density"
              Natmos = nHeps(el,0)*zweight(0)
              print*,"column densities:",Natmos,REAL(crust_Neps(el))
              if (ABS(crust_Neps(el))>0.5*Natmos) then
                print*,"*** too large timestep" 
                toomuch = .true. 
              endif  
            endif  
          endif  
        enddo  

        if (toomuch) then
          dt = dt/dt_increase
          nHeps(:,:) = nHold(:,:) 
          crust_Neps(:) = crust_Nold(:)
          reduced = .true.
          if (verbose>0) read'(A1)',char1
          goto 100
        endif  

        time = time + dt
        if (verbose>0) write(*,'(TL10," DIFFUSION:",I8,A,$)') it,CR
        if (time.ge.deltat) exit
        if (exhausted.or.reduced) exit

      enddo  

      deltat = time    ! return actually advanced timestep
      if (verbose>0) print'(" EXPLICIT DIFFUSION:",I8,"  time=",
     >               1pE11.4,"  Dt=",1pE11.4)',it,time0,time
      end


************************************************************************
      subroutine DIFFUSION_IMPLICIT(time0,deltat,isort,ipass,esort,
     >                              in_crust,limiting,reduced,spnam)
************************************************************************
      use GRID,ONLY: N=>Npoints,zz,zweight,d1l,d1m,d1r,d2l,d2m,d2r,
     >               dd1l,dd1m,dd1r,dt_diff_im,xlower,xupper
      use PARAMETERS,ONLY: bc_low,bc_high,verbose,
     >                     outflux,detailed_H
      use STRUCT,ONLY: Diff,nHtot,nHeps,crust_Neps,crust_Ncond
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use ELEMENTS,ONLY: NELEM,elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      use JEANS_ESCAPE,ONLY: EXTRA,Hfrac,jpern
      use EXCHANGE,ONLY: H
      implicit none
      real*8,intent(IN) :: time0
      real*8,intent(INOUT) :: deltat
      logical,intent(in) :: limiting(NELEM+EXTRA),in_crust(NELEM+EXTRA)
      integer,intent(in) :: isort,ipass,esort(NELEM+EXTRA)
      logical,intent(out) :: reduced
      character(len=3),intent(in) :: spnam(NELEM+EXTRA)
      integer :: i,j,it,e,el,Nstep,round
      real*8,dimension(0:N) :: xx,xnew,xHtot,nHmerk,rest
      real*8,dimension(N+1,N+1) :: BB
      real*8,dimension(N+1,N+1,NELEM+EXTRA) :: BBel
      real*8 :: influx(NELEM+EXTRA),Hinflux,dNcol,Natmos,NHatmos
      real*8 :: nHold(NELEM+EXTRA,-2:N),crust_Nold(NELEM)
      real*8 :: nD,time,dt,xl,xm,xr
      character :: CR = CHAR(13)
      character(len=1) :: char1
      logical :: exhausted,toomuch

      reduced = .false.
      exhausted = .false.

      Nstep = 20
      time  = 0.d0
      dt    = deltat/Nstep
      do e=1,isort
        el = esort(e)
        if (in_crust(el).and.limiting(e)) then   
          bc_low = 1     ! fixed concentration inner boundary condition
        else                              
          bc_low = 2     ! fixed flux inner boundary condition
        endif
        call INIT_DIFFUSION(el,N,bc_low,bc_high,dt,BB)  
        BBel(:,:,el) = BB(:,:)
      enddo  
      !NHatmos = 0.d0
      !do i=0,N
      !  NHatmos = NHatmos + nHeps(H,i)*zweight(i)
      !enddo
      !print*,"NH vorher:",NHatmos
      
      do it=1,Nstep

        nHold(:,:) = nHeps(:,:)
        crust_Nold(:) = crust_Neps(:)
 100    continue
        print*,"dt=",dt
        round = 1
        influx(:) = 0.d0
        Hinflux   = 0.d0
        xHtot(:)  = 0.d0
        if (MINVAL(nHeps)<0.d0) then
          print*,"*** negative nHeps it=",it
          stop
        endif  
        toomuch = .false.

        do e=1,isort

          el = esort(e)
          if ((round==1).and.(.not.limiting(e))) then
            round = 2 
            call INFLUXES(influx,isort,ipass,esort,limiting) 
          endif  
          if (el<=NELEM) then
            xx(0:N) = nHeps(el,0:N)/nHtot(0:N)                ! regular element
          else  
            xx(0:N) = nHeps(H,0:N)/nHtot(0:N)*Hfrac(el-NELEM) ! H,H2,H2O
          endif  
          BB(:,:) = BBel(:,:,el) 

          !------------------------------
          ! ***  boundary conditions  ***
          !------------------------------
          rest = xx
          nD = nHtot(0)*Diff(0)  
          if (in_crust(el).and.limiting(e)) then   
            rest(0) = xlower(el)                 ! fixed concentration
            bc_low = 1
          else               
            bc_low = 2
            rest(0) = -influx(el)/nD/dd1l(0)     ! fixed flux
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
          !if (verbose>0) then
          !  print'(A3,99(1pE10.3))',spnam(el),xx(0:3)
          !  print'(3x,99(1pE10.3))',xnew(0:3)
          !endif  
          if (verbose>0) then
            print'(A3,99(1pE10.3))',spnam(el),xx(N-3:N)
            print'(3x,99(1pE10.3))',xnew(N-3:N)
          endif  

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

          if (round==2) then
            print*,"deviation "//spnam(el),xx(1)/xnew(1)-1.d0
            print'(A3,"  influx=",2(1pE12.4))',
     >           elnam(el),influx(el)
            if (ABS(xnew(1)/xx(1)-1.d0)>0.01) then
              print*,'*** too large deviation'
              toomuch = .true.
            endif  
          endif

          if (el==H) then
            nHmerk(0:N) = nHtot(0:N)*xnew(0:N)
            print'("outflux <H>=",3(1pE15.7))',
     >           nHmerk(N),jpern(H),nHmerk(N)*jpern(H)
          else if (el<=NELEM) then
            nHeps(el,0:N) = nHtot(0:N)*xnew(0:N)
          else
            nHeps(el,0:N) = nHtot(0:N)*xnew(0:N)
            xHtot(0:N) = xHtot(0:N) + xnew(0:N)         ! add together new hydrogen
            Hinflux = Hinflux + influx(el)
          endif  
          if (el==NELEM+EXTRA) then
            if (detailed_H) then 
              !print*,nHtot(0:N)*xHtot(0:N)/nHmerk(0:N) 
              nHeps(H,0:N) = nHtot(0:N)*xHtot(0:N)
              influx(H) = Hinflux
              outflux = SUM(nHeps(NELEM+1:NELEM+EXTRA,N)
     >                     *jpern(NELEM+1:NELEM+EXTRA))
              print'("outflux Hat=",3(1pE15.7))',nHeps(NELEM+1,N),
     >             jpern(NELEM+1),nHeps(NELEM+1,N)*jpern(NELEM+1)
              print'("outflux H2 =",3(1pE15.7))',nHeps(NELEM+2,N),
     >             jpern(NELEM+2),nHeps(NELEM+2,N)*jpern(NELEM+2)
              print'("outflux H2O=",3(1pE15.7))',nHeps(NELEM+3,N),
     >             jpern(NELEM+3),nHeps(NELEM+3,N)*jpern(NELEM+3)
              print'("outflux <H>=",3(1pE15.7))',nHeps(H,N),outflux
            else  
              nHeps(H,0:N) = nHmerk(0:N)
            endif
          endif  

        enddo  
        
        !----------------------------------------------------
        ! ***  check whether influx would substantially   *** 
        ! ***  increase total atmospheric column density  ***
        !----------------------------------------------------
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
              toomuch = .true.
            endif  
          endif
        enddo     

        do e=1,NELM
          if (e==iel) cycle
          el = elnum(e)
          !----------------------------------------
          ! ***  update crust column densities  ***
          !----------------------------------------
          if (in_crust(el)) then 
            if (verbose>0) then
              print*,elnam(el),REAL(crust_Neps(el)),influx(el)*dt
            endif  
            crust_Neps(el) = crust_Neps(el) - influx(el)*dt
            if (crust_Neps(el)<0.Q0) then
              exhausted = .true. 
              reduced = .true.
              print*,elnam(el)," *** negative crust column density"
              Natmos = nHeps(el,0)*zweight(0)
              print*,"column densities:",Natmos,REAL(crust_Neps(el))
              if (ABS(crust_Neps(el))>0.5*Natmos) then
                print*,"*** too large timestep" 
                toomuch = .true. 
              endif  
            endif  
          endif  
        enddo  

        if (toomuch) then
          dt = dt/2.0
          reduced = .true.
          nHeps(:,:) = nHold(:,:) 
          crust_Neps(:) = crust_Nold(:)
          do e=1,isort
            el = esort(e)
            if (in_crust(el).and.limiting(e)) then   
              bc_low = 1
            else                              
              bc_low = 2
            endif
            call INIT_DIFFUSION(el,N,bc_low,bc_high,dt,BB)  
            BBel(:,:,el) = BB(:,:)
          enddo  
          if (verbose>0) read'(A1)',char1
          goto 100
       endif  

        time = time + dt
        if (verbose>0) write(*,'(TL10," DIFFUSION:",I8,A,$)') it,CR
        if (exhausted.or.reduced) exit

      enddo

      !NHatmos = 0.d0
      !do i=0,N
      !  NHatmos = NHatmos + nHeps(H,i)*zweight(i)
      !enddo
      !print*,"NH nachher:",NHatmos

      deltat = time    ! return actually advanced timestep
      if (verbose>0) print'(" IMPLICIT DIFFUSION:",I8,"  time=",
     >               1pE11.4,"  Dt=",1pE11.4)',it,time0,time
      end


************************************************************************
      SUBROUTINE INFLUXES(influx,isort,ipass,esort,limiting) 
************************************************************************
      use PARAMETERS,ONLY: verbose
      use STRUCT,ONLY: crust_Neps,crust_Ncond
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_nel,dust_el,dust_nu
      use ELEMENTS,ONLY: NELEM,elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      use JEANS_ESCAPE,ONLY: EXTRA
      implicit none
      real*8,intent(inout) :: influx(NELEM+EXTRA)
      logical,intent(in) :: limiting(NELEM+EXTRA)
      integer,intent(in) :: isort,ipass,esort(NELEM+EXTRA)
      integer :: i,ii,jj,e,el,el2,i1,i2,s,Neq1,Neq2
      real*8 :: AA(NELM,NELM),sol(NELM),rhs(NELM)
      real*8 :: dNcond(NDUST),check(NELEM+EXTRA),qual

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
        if (el>NELEM) cycle
        ii = ii+1
        rhs(ii) = influx(el)
        AA(ii,:) = 0.d0
        jj = 0
        do s=1,NDUST
          if (crust_Ncond(s)<=0.Q0) cycle
          jj = jj+1
          if (e==i1.and.verbose>1) print*,jj,dust_nam(s)
          Neq2 = jj
          do i=1,dust_nel(s)
            el2 = dust_el(s,i)
            if (el==el2) AA(ii,jj) = dust_nu(s,i)
          enddo
        enddo
        if (verbose>1) print*,ii,elnam(el)
      enddo  
      if (verbose>1) then
        print*,Neq1,Neq2
        do i=1,Neq1
          print'(99(I3))',int(AA(i,1:Neq2))
        enddo  
      endif  
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
        if (el>NELEM) cycle
        do s=1,NDUST
          if (crust_Ncond(s)<=0.Q0) cycle
          do i=1,dust_nel(s)
            el2 = dust_el(s,i)
            if (el==el2) then
              influx(el) = influx(el)+dust_nu(s,i)*dNcond(s)
            endif  
          enddo
        enddo
        if (verbose>0) then
          if (limiting(e)) then
            qual = ABS(influx(el)/check(el)-1.d0) 
            print'(A3,"  influx=",1pE12.4,"  check=",1pE9.2)',
     >                   elnam(el),influx(el),qual
          else
            print'(A3,"  influx=",1pE12.4)',elnam(el),influx(el)
          endif  
        endif  
      enddo  

      end
