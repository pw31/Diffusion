**********************************************************************
      MODULE DATABASE
**********************************************************************
      use ELEMENTS,ONLY: NELEM
      use DUST_DATA,ONLY: NDUSTmax,dust_nam
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,parameter :: DMAX = 2*10**5
      integer :: NDAT=0,NLAST=0,NMODI=0,NPICK1=0,NPICK2=0
      TYPE ENTRY
        real*8 :: ln
        real*8 :: lT
        real*8 :: sumeps
        real(kind=qp) :: eps(NELEM)
        real(kind=qp) :: ddust(NDUSTmax)
        integer :: Nsolve
        integer :: indep(NELEM)
      END TYPE ENTRY
      TYPE(ENTRY) :: dbase(DMAX)
      end MODULE DATABASE

**********************************************************************
      SUBROUTINE SAVE_DBASE
**********************************************************************
      use ELEMENTS,ONLY: NELEM
      use DUST_DATA,ONLY: NDUST,dust_nam
      use DATABASE,ONLY: NDAT,NLAST,dbase
      implicit none
      integer :: i
      character(len=80) :: filename="database.dat"
      if (NLAST==0) then
        open(unit=11,file=filename,form='unformatted',status='replace')
        write(11) NELEM,NDUST
        write(11) dust_nam
        do i=1,NDAT
          write(11) dbase(i)%ln 
          write(11) dbase(i)%lT
          write(11) dbase(i)%sumeps
          write(11) dbase(i)%eps
          write(11) dbase(i)%ddust(1:NDUST)
          write(11) dbase(i)%Nsolve
          write(11) dbase(i)%indep(1:dbase(i)%Nsolve)
        enddo 
        close(11)
      else if (NDAT>NLAST) then 
        open(unit=11,file=filename,form='unformatted',position='append')
        do i=NLAST+1,NDAT
          write(11) dbase(i)%ln 
          write(11) dbase(i)%lT
          write(11) dbase(i)%sumeps
          write(11) dbase(i)%eps
          write(11) dbase(i)%ddust(1:NDUST)
          write(11) dbase(i)%Nsolve
          write(11) dbase(i)%indep(1:dbase(i)%Nsolve)
        enddo 
        close(11)
      endif  
      NLAST = NDAT
      end

**********************************************************************
      SUBROUTINE LOAD_DBASE
**********************************************************************
      use PARAMETERS,ONLY: verbose
      use ELEMENTS,ONLY: NELEM
      use DUST_DATA,ONLY: NDUST,dust_nam
      use DATABASE,ONLY: qp,NDAT,NLAST,dbase
      implicit none
      integer :: i,NELEM_read,NDUST_read
      logical :: ex
      character(len=20) :: dust_nam_read(NDUST)
      character(len=80) :: filename="database.dat"

      NDAT = 0
      NLAST = 0
      inquire(file=filename,exist=ex)
      if (.not.ex) goto 200
      open(unit=11,file=filename,form="unformatted",status="old")
      read(11) NELEM_read,NDUST_read
      if (NELEM_read.ne.NELEM) goto 200
      if (NDUST_read.ne.NDUST) goto 200
      read(11) dust_nam_read
      do i=1,NDUST
        if (dust_nam(i).ne.dust_nam_read(i)) goto 200
      enddo
      do i=1,999999
        read(11,end=100) dbase(i)%ln 
        read(11) dbase(i)%lT
        read(11) dbase(i)%sumeps
        read(11) dbase(i)%eps
        read(11) dbase(i)%ddust(1:NDUST)
        read(11) dbase(i)%Nsolve
        read(11) dbase(i)%indep(1:dbase(i)%Nsolve)
        NDAT = NDAT+1
      enddo 
 100  close(11)
      if (verbose>0) print*,"... having read ",NDAT," datasets." 
      NLAST = NDAT
      return
 200  close(11)
      if (verbose>0) print*,"... no / unsuitable database."
      end

**********************************************************************
      SUBROUTINE PUT_DATA(nH,T,eps,ddust,qbest,ibest,active,
     >                    Nsolve,indep)
**********************************************************************
      use PARAMETERS,ONLY: verbose
      use ELEMENTS,ONLY: NELEM,eps0
      use DUST_DATA,ONLY: NDUST
      use DATABASE,ONLY: qp,NDAT,NLAST,NMODI,DMAX,dbase
      implicit none
      real*8,intent(in) :: nH,T,qbest
      integer,intent(in) :: ibest,Nsolve,indep(Nsolve)
      real(kind=qp),intent(in) :: eps(NELEM),ddust(NDUST)
      real*8 :: sumeps
      logical,intent(in) :: active(0:NDUST)
      integer :: i,j
      
      sumeps = 0.0
      do i=1,NELEM
        sumeps = sumeps+eps0(i)
      enddo  
      if (qbest<1.d-8) then
        return 
      else if (ibest>0.and.qbest<1.d-4) then
        i = ibest
        if (verbose>0) then
          write(*,'(" ... replacing database entry (",I6,
     >          ") nH,T,sumeps=",3(1pE15.7))') i,nH,T,sumeps
        endif  
      else  
        NDAT = NDAT+1
        i = NDAT
        if (verbose>0) then
          write(*,'(" ... adding database entry (",I6,
     >          ") nH,T,sumeps=",3(1pE15.7))') i,nH,T,sumeps
        endif  
        if (NDAT>DMAX) then
          print*,"*** NDAT>DMAX in PUT_DATA",NDAT,DMAX
          stop
        endif  
      endif  
      dbase(i)%ln     = LOG(nH)
      dbase(i)%lT     = LOG(T)
      dbase(i)%sumeps = sumeps
      dbase(i)%eps    = eps
      do j=1,NDUST
        dbase(i)%ddust(j) = ddust(j)
        if (.not.active(j)) dbase(i)%ddust(j)=0.Q0
      enddo
      dbase(i)%Nsolve = Nsolve
      dbase(i)%indep(1:Nsolve) = indep(1:Nsolve)
      NMODI = i
      if (NDAT>NLAST+10) then
        call SAVE_DBASE
        print*,"... saved ",NDAT," datasets."
      endif  
      end


**********************************************************************
      subroutine GET_DATA(nH,T,eps,ddust,qbest,ibest,active,verbose)
**********************************************************************
      use ELEMENTS,ONLY: NELEM,NEPS,eps0,elnam,elnr
      use dust_data,ONLY: NDUST,dust_nel,dust_nu,dust_el,dust_nam
      use DATABASE,ONLY: qp,NDAT,NMODI,NPICK1,NPICK2,DMAX,dbase
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      implicit none
      real*8,intent(in) :: nH,T
      real*8,intent(out) :: qbest
      integer,intent(out) :: ibest
      real(kind=qp),intent(inout) :: eps(NELEM),ddust(NDUST)
      logical,intent(out) :: active(0:NDUST)
      integer,intent(in) :: verbose
      real*8 :: ln,lT,seps,lnread,lTread,sepsread,qgood,qual,pot
      real*8 :: rsort(NEPS),sort
      real(kind=qp) :: check(NELEM),error,errmax,neweps,dd,dmin
      real(kind=qp) :: stoich(NEPS,NEPS),xx(NEPS),rest(NEPS)
      integer :: i,j,k,it,el,elworst,Ncond,Nact,sp,iloop,imin
      integer :: isort(NEPS),dact(NDUST),OK,ecount(NELEM)
      integer :: Nsolve,indep(NELEM)
      character(len=1) :: char
      logical :: found,eact(NELEM),epure(NELEM),eflag(NELEM)
      logical :: nan,IS_NAN
      logical,save :: firstCall=.true.
      
      if (firstCall) then
        call LOAD_DBASE 
        firstCall = .false.
      endif  

      seps = 0.0
      do i=1,NELEM
        seps = seps + eps0(i)
      enddo  
      if (verbose>0) then
        print'("looking for nH,T,sumeps=",3(1pE13.5)," ...")',nH,T,seps
      endif  
      ln = LOG(nH)
      lT = LOG(T) 
      qbest  = 9.d+99
      ibest  = 0
      pot    = -0.03
      qgood  = 1.E-5
      !--- try last entry modified first ---
      if (NMODI>0) then
        i=NMODI
        lnread   = dbase(i)%ln 
        lTread   = dbase(i)%lT
        sepsread = dbase(i)%sumeps
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
     >       + 1.E-5*ABS(sepsread-seps)
        qbest = qual
        ibest = i
        if (qbest<qgood) goto 100
      endif  
      !--- try around entry picked last time ---  
      do i=MAX(1,NPICK1-1),MIN(NDAT,NPICK1+1)
        lnread   = dbase(i)%ln 
        lTread   = dbase(i)%lT
        sepsread = dbase(i)%sumeps
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
     >       + 1.E-5*ABS(sepsread-seps)
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<qgood) goto 100
        endif  
      enddo
      do i=MAX(1,NPICK2-1),MIN(NDAT,NPICK2+1)
        lnread   = dbase(i)%ln 
        lTread   = dbase(i)%lT
        sepsread = dbase(i)%sumeps
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
     >       + 1.E-5*ABS(sepsread-seps)
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<qgood) goto 100
        endif  
      enddo
      if (verbose>0) then
        write(*,*) "entering full search ..."
      endif  
      !--- check them all ---  
      do i=NDAT,1,-1
        lnread   = dbase(i)%ln 
        lTread   = dbase(i)%lT
        sepsread = dbase(i)%sumeps
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
     >       + 1.E-5*ABS(sepsread-seps)
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<qgood) goto 100
        endif  
      enddo
 100  active = .false.
      if (ibest==0) return

      if (verbose>0) then
        print'(" ... found best dataset (",I6,")  nH,T,sumeps,qual=",
     >       4(1pE13.5))',ibest,EXP(dbase(ibest)%ln),
     >       EXP(dbase(ibest)%lT),dbase(ibest)%sumeps,qbest
      endif  
      NPICK2 = NPICK1
      NPICK1 = ibest
      eps    = dbase(ibest)%eps
      ddust  = dbase(ibest)%ddust(1:NDUST)
      Nsolve = dbase(ibest)%Nsolve
      indep(1:Nsolve) = dbase(ibest)%indep(1:Nsolve)
      
      do iloop=1,5
        !--- 1. determine active condensates and involved elements ---
        Ncond = 0
        eact = .false.
        epure = .false.
        ecount = 0
        do i=1,NDUST
          if (ddust(i)>0.Q0) then
            active(i)=.true.
            Ncond = Ncond+1
            dact(Ncond) = i
            do j=1,dust_nel(i)
              el = dust_el(i,j) 
              eact(el) = .true.
              ecount(el) = ecount(el)+1
              if (dust_nel(i)==1) epure(el)=.true.
            enddo
          endif  
        enddo
        !--- 2. check element conservation ---
        check = eps
        do i=1,NDUST
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            check(el) = check(el) + ddust(i)*dust_nu(i,j)    
          enddo
        enddo
        errmax = -1.Q0
        do i=1,NELM
          if (i==iel) cycle
          el = elnum(i)
          error = ABS(1.Q0-check(el)/eps0(el))
          if (error.gt.errmax) then
            errmax = error
            elworst = el
          endif   
        enddo
        if (verbose>1) print*,"worst element conservation "
     >                      //elnam(elworst),errmax
        if (errmax<1.Q-25) return ! perfect fit - nothing to do
        !--- 3. sort elements according to gas abundance ---
        rsort = 9.d+99
        isort = 0
        Nact  = 0
        do i=1,NEPS
          el = elnr(i)
          if (.not.eact(el)) cycle
          sort = eps(el)
          if (epure(el).and.ecount(el)==1) sort=0.0
          do j=NEPS+1,2,-1
            if (sort>rsort(j-1)) exit
          enddo  
          isort(j+1:NEPS) = isort(j:NEPS-1)
          rsort(j+1:NEPS) = rsort(j:NEPS-1)
          isort(j) = el
          rsort(j) = sort
          Nact = Nact+1
        enddo
        do j=1,Nact
          el = isort(j)
          if (.not.eact(el)) cycle
          if (verbose>1) print'(A3,2(1pE12.5))',
     >                   elnam(el),eps0(el),eps(el)
        enddo
        if (verbose>1) then
          do i=1,Ncond
            print*,i,dust_nam(dact(i)) 
          enddo  
          print'("number of active condensates and elements =",2(I3))',
     >       Ncond,Nact
        endif  
        !--- take sorting from database if possible ---
        if (verbose>1) then
          print*,Nsolve,Ncond
          print*,elnam(isort(1:Ncond))
          print*,elnam(indep(1:Nsolve))
        endif  
        if (Nsolve==Ncond) then
          isort(1:Nsolve) = indep(1:Nsolve)
          eflag(:) = .true.
          eflag(isort(1:Nsolve)) = .false.
          j=Nsolve
          do i=1,NELM
            if (i==iel) cycle
            el = elnum(i)
            if (.not.eact(el)) cycle
            if (eflag(el)) then
              j = j + 1 
              isort(j) = el
              if (verbose>1) then
                print*,j,elnam(el)//" non-limiting element"
              endif  
            endif
          enddo
        endif  
        !--- 4. find fitting linear combination of solids ---
        stoich = 0.Q0
        rest = 0.Q0
        do j=1,Ncond
          el = isort(j)
          rest(j) = eps0(el)-eps(el)
          do i=1,Ncond 
            sp = dact(i)
            do k=1,dust_nel(sp)
              if (dust_el(sp,k)==el) then
                stoich(j,i) = dust_nu(sp,k) 
              endif
            enddo
          enddo
        enddo  
        if (verbose>1) then
          do i=1,Ncond
            print'(99(I3))',int(stoich(i,1:Ncond)) 
          enddo  
        endif  
        !-------------------------------------------
        call GAUSS16( NEPS, Ncond, stoich, xx, rest)
        !-------------------------------------------
        OK = 0
        imin = 0
        dmin = 9.e+99
        do i=1,Ncond
          sp = dact(i)
          dd = 0.0
          do k=1,dust_nel(sp)
            el = dust_el(sp,k)
            dd = MAX(dd,dust_nu(sp,k)*ddust(sp)/(eps0(el)-eps(el)))
          enddo  
          if (verbose>1) print'(A15,1pE9.2,1pE16.8," ->",1pE16.8)',
     >                       dust_nam(sp),dd,ddust(sp),xx(i)
          if (IS_NAN(DBLE(xx(i)))) then
            OK = -4
            exit
          endif  
          if (xx(i)<0.Q0) then
            OK = -1
            if (dd<dmin) then
              dmin = dd 
              imin = sp
            endif  
          endif  
        enddo
        if (OK==0) then
          do i=1,Ncond
            sp = dact(i)
            ddust(sp) = xx(i)  
          enddo  
          exit
        endif
        if (OK==-4) exit
        active(imin) = .false.    
        ddust(imin) = 0.Q0
        print*,"switching off "//trim(dust_nam(imin))//" ..."
      enddo 
      if (OK.ne.0) goto 200
      !--- 5. correct non-involved and dependent element abundances ---
      do el=1,NELEM
        if (.not.eact(el)) then
          eps(el)=eps0(el)
          !print*,"non-involved ",elnam(el),eps0(el)
        endif  
      enddo
      do j=Ncond+1,Nact
        el = isort(j)
        neweps = eps0(el)
        do i=1,Ncond 
          sp = dact(i)
          do k=1,dust_nel(sp)
            if (dust_el(sp,k)==el) then
              neweps = neweps - dust_nu(sp,k)*ddust(sp) 
            endif
          enddo
        enddo
        if (verbose>1) print'("dependent",A3,1pE16.8," ->",1pE16.8)',
     >                      elnam(el),eps(el),neweps
        eps(el) = neweps
      enddo  
      !--- 6. check element conservation ---
      check = eps
      do i=1,NDUST
        do j=1,dust_nel(i)
          el = dust_el(i,j)
          check(el) = check(el) + ddust(i)*dust_nu(i,j)    
        enddo
      enddo
      errmax = -1.Q0
      do i=1,NELM
        if (i==iel) cycle
        el = elnum(i)
        error = ABS(1.Q0-check(el)/eps0(el))
        if (verbose>1) print'(A15,2(1pE12.5))',elnam(el),eps(el),error
        if (error.gt.errmax) then
          errmax = error
          elworst = el
        endif   
        if (eps(el)<0.Q0) OK=-2
      enddo
      !print*,"worst element conservation "//elnam(elworst),errmax
      if (errmax>1.Q-20) OK=-3

 200  continue
      if (OK==-1) then
        print*,"*** negative dust abundances in database.f"
      else if (OK==-2) then
        print*,"*** negative element abundances in database.f"
      else if (OK==-3) then
        print*,"*** element conservation error in database.f"
      else if (OK==-4) then
        print*,"*** NaN in database.f"
      endif
      if (OK<0) then
        ddust(j) = 0.Q0
        active(j) = .false.
        qbest = 9.d+99
      endif  
      if (verbose>2) read(*,'(A1)') char
          
      end
