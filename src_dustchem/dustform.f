************************************************************************
      SUBROUTINE DUSTFORM(time,deltat,verbose)
************************************************************************
      use GRID,ONLY: Npoints
      use STRUCT,ONLY: nHtot,Temp,nHeps,rhoLj,rhoL3
      use ELEMENTS,ONLY: NEPS,elnam,elcode
      use CHEMISTRY,ONLY: NMOLE
      use DUST_DATA,ONLY: NDUST,dust_Vol,dust_nam,
     >                    dust_nel,dust_el,dust_nu
      use NUCLEATION,ONLY: NNUC
      use EXCHANGE,ONLY: ipoint,Jst,chi,nmol
      implicit none
      integer,intent(in) :: verbose
      real*8,intent(in) :: time,deltat
      real*8,allocatable :: yy(:),FF(:)
      real*8 :: tt,dt,corr,y0,y1,y2,y3,bb,cc,epsd,nHeldust,stoich
      integer :: i,j,ip,el,NN,ifail
      logical :: has_dust,has_rate,evap
      character(len=1) :: char1
      character :: CR=CHAR(13)

      if (verbose>1) then
        print*
        print*,"entering DUSTFORM ..."
        print*,"====================="
      endif  
      NN = 4+NDUST+NEPS

!$omp parallel 
!$omp& default(none)
!$omp& shared(NN,NNUC,NDUST,NPOINTS,NEPS,NMOLE,CR,dust_Vol,verbose)
!$omp& shared(rhoLj,rhoL3,nHtot,nHeps,Temp,time,deltat,elcode,elnam)
!$omp& shared(dust_nam,dust_nel,dust_el,dust_nu)
!$omp& private(i,j,ip,yy,FF,has_dust,has_rate,evap,char1)
!$omp& private(tt,dt,corr,y0,y1,y2,y3,bb,cc,epsd,ifail)
!$omp& private(el,nHeldust,stoich)
      allocate(yy(NN),FF(NN))
      if (.not.allocated(nmol)) then
        allocate(nmol(NMOLE),Jst(NNUC),chi(NDUST))
      endif

!$omp do schedule(dynamic,1)
      do ip=2,Npoints-1        ! ip=1 and Npoints will be overrules by
                               ! diffusion boundary conditions anyway
        ipoint = ip 
        yy(1:4) = rhoLj(0:3,ipoint)            ! rho*Lj    (j=0,1,2,3)
        yy(5:4+NDUST) = rhoL3(1:NDUST,ipoint)  ! rho*L3^s  (s=1...NDUST)
        yy(5+NDUST:NN) = nHeps(1:NEPS,ipoint)  ! element abundances
        epsd = rhoLj(3,ipoint)/dust_Vol(1)/nHtot(ipoint)

        tt = 0.d0
        call DUST_RHS(NN,tt,yy,FF,.false.,ifail) 
        has_dust = ((MAXVAL(chi)>0.0.and.epsd>1.E-15).or.(epsd>1.E-10))
        has_rate = (SUM(Jst)>0.0)
        evap = .false.

        !if (yy(1)>0.d0.and.FF(4)<0.d0) then
        !  !--- smooth it a bit ---   
        !  y0 = LOG(yy(1))
        !  y1 = LOG(yy(2))
        !  y2 = LOG(yy(3))
        !  y3 = LOG(yy(4))
        !  cc = 0.25*(y0-y1-y2+y3)
        !  bb = (y3-y0)/3.0-3.0*cc
        !  !y1 = y0 + bb + cc
        !  !y2 = y0 + 2.0*bb + 4.0*cc
        !  !y3 = y0 + 3.0*bb + 9.0*cc
        !  if (verbose>1) print*,"--- smooth Lj ---"
        !  if (verbose>1) print'(4(1pE12.4))',yy(1:4)
        !  yy(2) = SQRT(yy(2) * EXP(y0 + bb + cc))
        !  yy(3) = SQRT(yy(3) * EXP(y0 + 2.0*bb + 4.0*cc))
        !  if (verbose>1) print'(4(1pE12.4))',yy(1:4)
        !  call DUST_RHS(NN,tt,yy,FF,.false.,ifail)
        !endif  

!$omp critical(output)
        if (verbose>1) then
          print*,"----------------------------------------------------"
          !print'(" yy=",99(1pE9.2))',yy
          !print'(" FF=",99(1pE9.2))',FF
          print'(I4," n<H>=",1pE9.2,"  T=",0pF8.2,
     >         "  has_dust=",L1,"  has_rate=",L1,"  epsd=",1pE9.2)',
     >         ipoint,nHtot(ipoint),Temp(ipoint),has_dust,has_rate,epsd
        endif
        if (verbose==1) write(*,'(TL10," DUSTFORM: ",I8,A,$)') ipoint,CR
!$omp end critical(output)

        if (has_dust.or.has_rate) then 
          tt = 0.d0
          dt = 1.e+99
          do j=1,NN
            if (yy(j)>0.d0.and.ABS(FF(j))>1.d-99) then
              dt = MIN(dt,0.01*yy(j)/ABS(FF(j)))
            endif  
          enddo  
          if (verbose>1) print'(" dt=",1pE9.2)',dt
          if (dt>1.E+30*deltat) goto 100
          if (verbose>1) call DUST_RHS(NN,tt,yy,FF,.true.,ifail)

          if (yy(1)==0.d0) then
            !--- two tiny explicit steps to kickstart dust formation --- 
            yy = yy + 1.d-10*deltat*FF
            tt = tt + 1.d-10*deltat
            call DUST_RHS(NN,tt,yy,FF,.false.,ifail)
            yy = yy + 1.d-10*deltat*FF
            tt = tt + 1.d-10*deltat
          endif  
          call DUST_TSTEP(NN,yy,FF,tt,dt,deltat,verbose,evap)
          if (verbose>1) call DUST_RHS(NN,tt,yy,FF,.true.,ifail) 
        else
          evap = .true. 
        endif  

 100    continue
!$omp critical(update)
        nHeps(1:NEPS,ipoint) = yy(5+NDUST:NN)        ! element abundances
        if (evap) then
          if (has_dust.and.verbose>1) print*," all dust evaporated"
          do i=1,NDUST
            nHeldust = rhoL3(i,ipoint)/dust_Vol(i)
            !print*,dust_nam(i),nHeldust/nHtot(ipoint)
            do j=1,dust_nel(i)            
              el = dust_el(i,j)
              stoich = dust_nu(i,j)
              !print*,elnam(el),stoich
              el = elcode(el)
              nHeps(el,ipoint) = nHeps(el,ipoint) + nHeldust*stoich
            enddo
          enddo  
          rhoLj(0:3,ipoint) = 0.d0                   ! rho*Lj    (j=0,1,2,3)
          rhoL3(1:NDUST,ipoint) = 0.d0               ! rho*L3^s  (s=1...NDUST)
        else  
          corr = 0.d0 
          do i=5,4+NDUST
            corr = corr + yy(i)
          enddo
          if (corr.gt.0.d0) then
            corr = yy(4)/corr
          else
            corr = 1.d0
          endif  
          !print*,"corr=",corr
          rhoLj(0:3,ipoint) = yy(1:4)                ! rho*Lj    (j=0,1,2,3)
          rhoL3(1:NDUST,ipoint) = yy(5:4+NDUST)*corr ! rho*L3^s  (s=1...NDUST)
        endif  
        if ((has_dust.or.has_rate).and.verbose>1) read(*,'(A1)') char1
!$omp end critical(update)

      enddo
!$omp end parallel      

      if (verbose==1) print'(" DUSTFORM: "I8,"  time=",1pE11.4,
     >                "  Dt=",1pE11.4)',Npoints,time,deltat

      end
