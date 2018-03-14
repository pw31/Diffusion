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
        real(kind=qp) :: eps(NELEM)
        real(kind=qp) :: ddust(NDUSTmax)
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
          write(11) dbase(i)%eps
          write(11) dbase(i)%ddust(1:NDUST)
        enddo 
        close(11)
      else if (NDAT>NLAST) then 
        open(unit=11,file=filename,form='unformatted',position='append')
        do i=NLAST+1,NDAT
          write(11) dbase(i)%ln 
          write(11) dbase(i)%lT
          write(11) dbase(i)%eps
          write(11) dbase(i)%ddust(1:NDUST)
        enddo 
        close(11)
      endif  
      NLAST = NDAT
      end

**********************************************************************
      SUBROUTINE LOAD_DBASE
**********************************************************************
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
        read(11) dbase(i)%eps
        read(11) dbase(i)%ddust(1:NDUST)
        NDAT = NDAT+1
        !print*,i,EXP(dbase(i)%ln),EXP(dbase(i)%lT)
      enddo 
 100  close(11)
      print*,"... having read ",NDAT," datasets." 
      NLAST = NDAT
      return
 200  close(11)
      print*,"... no / unsuitable database."
      end

**********************************************************************
      SUBROUTINE PUT_DATA(nH,T,eps,ddust,qbest,ibest,active)
**********************************************************************
      use ELEMENTS,ONLY: NELEM
      use DUST_DATA,ONLY: NDUST
      use DATABASE,ONLY: qp,NDAT,NLAST,NMODI,DMAX,dbase
      implicit none
      real*8,intent(in) :: nH,T,qbest
      integer,intent(in) :: ibest
      real(kind=qp),intent(in) :: eps(NELEM),ddust(NDUST)
      logical,intent(in) :: active(0:NDUST)
      integer :: i,j
      
      if (qbest<1.d-8) then
        return 
      else if (ibest>0.and.qbest<1.d-3) then
        i = ibest
        write(*,'(" ... replacing database entry (",I6,
     >          ") nH,T=",2(1pE15.7))') i,nH,T
      else  
        NDAT = NDAT+1
        i = NDAT
        write(*,'(" ... adding database entry (",I6,
     >          ") nH,T=",2(1pE15.7))') i,nH,T
        if (NDAT>DMAX) then
          print*,"*** NDAT>DMAX in PUT_DATA",NDAT,DMAX
          stop
        endif  
      endif  
      dbase(i)%ln = LOG(nH)
      dbase(i)%lT = LOG(T)
      dbase(i)%eps = eps
      do j=1,NDUST
        dbase(i)%ddust(j) = ddust(j)
        if (.not.active(j)) dbase(i)%ddust(j)=0.Q0
      enddo
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
      real*8 :: ln,lT,lnread,lTread,qual,pot,rsort(NEPS)
      real(kind=qp) :: check(NELEM),error,errmax,neweps
      real(kind=qp) :: stoich(NEPS,NEPS),xx(NEPS),rest(NEPS)
      integer :: i,j,k,it,el,elworst,Ncond,Nact,sp
      integer :: isort(NEPS),dact(NDUST),OK
      character(len=1) :: char
      logical :: found,eact(NELEM)
      logical,save :: firstCall=.true.
      
      if (firstCall) then
        call LOAD_DBASE 
        firstCall = .false.
      endif  

      write(*,'("looking for nH,T=",2(1pE13.5)," ...")') nH,T
      ln = LOG(nH)
      lT = LOG(T) 
      qbest  = 9.d+99
      ibest  = 0
      pot    = -0.03
      !--- try last entry modified first ---
      if (NMODI>0) then
        i=NMODI
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
        qbest = qual
        ibest = i
        if (qbest<1.d-3) goto 100
      endif  
      !--- try around entry picked last time ---  
      do i=MAX(1,NPICK1-1),MIN(NDAT,NPICK1+1)
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<1.d-3) goto 100
        endif  
      enddo
      do i=MAX(1,NPICK2-1),MIN(NDAT,NPICK2+1)
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<1.d-3) goto 100
        endif  
      enddo
      write(*,*) "entering full search ..."
      !--- check them all ---  
      do i=NDAT,1,-1
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        qual = 0.05*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<1.d-3) goto 100
        endif  
      enddo
 100  active = .false.
      if (ibest>0) then
        write(*,'(" ... found best dataset (",I6,
     >          ")  nH,T,qual=",3(1pE13.5))')
     >     ibest,EXP(dbase(ibest)%ln),EXP(dbase(ibest)%lT),qbest
        NPICK2 = NPICK1
        NPICK1 = ibest
        eps    = dbase(ibest)%eps
        ddust  = dbase(ibest)%ddust(1:NDUST)

        !--- 1. determine active condensates and involved elements ---
        Ncond = 0
        eact = .false.        
        do i=1,NDUST
          if (ddust(i)>0.Q0) then
            active(i)=.true.
            Ncond = Ncond+1
            dact(Ncond) = i
            do j=1,dust_nel(i)
              eact(dust_el(i,j))=.true.
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
        print*,"worst element conservation "//elnam(elworst),errmax
        if (errmax<1.Q-25) return          ! perfect fit - nothing to do
        !--- 3. sort elements according to gas abundance ---
        rsort = 9.d+99
        isort = 0
        Nact  = 0
        do i=1,NEPS
          el = elnr(i)
          if (.not.eact(el)) cycle
          do j=NEPS+1,2,-1
            if (eps(el)>rsort(j-1)) exit
          enddo  
          isort(j+1:NEPS) = isort(j:NEPS-1)
          rsort(j+1:NEPS) = rsort(j:NEPS-1)
          isort(j) = el
          rsort(j) = eps(el)
          Nact = Nact+1
        enddo
        do j=1,Nact
          el = isort(j)
          if (.not.eact(el)) cycle
          print'(A3,2(1pE12.5))',elnam(el),eps0(el),eps(el)
        enddo
        do i=1,Ncond
          print*,i,dust_nam(dact(i)) 
        enddo  
        print'("number of active condensates and elements =",2(I3))',
     >       Ncond,Nact
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
        do i=1,Ncond
          print'(99(I3))',int(stoich(i,1:Ncond)) 
        enddo  
        !-------------------------------------------
        call GAUSS16( NEPS, Ncond, stoich, xx, rest)
        !-------------------------------------------
        OK = 0
        do i=1,Ncond
          sp = dact(i)
          print'(A15,1pE16.8," ->",1pE16.8)',
     >          dust_nam(sp),ddust(sp),xx(i)
          ddust(sp) = xx(i)
          if (xx(i)<0.Q0) OK=-1
        enddo
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
          print'("dependent",A3,1pE16.8," ->",1pE16.8)',
     >         elnam(el),eps(el),neweps
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
          print'(A15,2(1pE12.5))',elnam(el),eps(el),error
          if (error.gt.errmax) then
            errmax = error
            elworst = el
          endif   
          if (eps(el)<0.Q0) OK=-2
        enddo
        !print*,"worst element conservation "//elnam(elworst),errmax
        if (errmax>1.Q-20) OK=-3

        if (OK==-1) then
          print*,"*** negative dust abundances in database.f"
        else if (OK==-2) then
          print*,"*** negative element abundances in database.f"
        else if (OK==-3) then
          print*,"*** element conservation error in database.f"
        endif
        if (OK<0) then
          ddust(j) = 0.Q0
          active(j) = .false.
          qbest = 9.d+99
          return
        endif  
        read(*,'(A1)') char
          
      endif
  
      end
