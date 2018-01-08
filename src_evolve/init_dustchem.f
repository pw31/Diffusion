**********************************************************************
      SUBROUTINE INIT_DUSTCHEM
**********************************************************************
      use PARAMETERS,ONLY: dustchem_file
      use NATURE,ONLY: amu
      use ELEMENTS,ONLY: NELEM,NEPS,mass,elnr,elcode,elnam
      use CHEMISTRY,ONLY: NMOLE,NELM,catm,cmol
      use DUST_DATA,ONLY: NDUSTmax,NDUST,Tmelt,Tcorr,
     &                    dust_nam,dust_rho,dust_vol,dust_mass,
     &                    dust_nel,dust_nu,dust_el,fit,cfit
      use NUCLEATION,ONLY: NNUC,nuc,nuc_nam
      use SPECIES,ONLY: NSPECIES,spnam,spnr,spzahl,spelem,spstoich,
     &                  keysp,spmass
      use REACTIONS,ONLY: NREAC,neduct,nprod,reac_sp,reac_nu,stick
      implicit none
      integer :: i,imax,j,k,el,j1,j2,nn,mm,kk,sp,stindex
      integer :: elem_cons(NELEM)
      real*8 :: dmass,prec(NDUSTmax)
      character(len=10000) :: allcond
      character(len=200):: zeile,lastzeile
      character(len=100) :: trivial(NDUSTmax),tmp
      character(len=20) :: name20
      character(len=2)  :: name2
      logical :: is_atom,found,allfound

      write(*,*) 
      write(*,*) "reading "//trim(dustchem_file)//" ..."
      write(*,*) "========================"
      trivial(:)=' '

      open(12, file='data/'//trim(dustchem_file), status='old')
 
      write(*,*) '--- dust species ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) imax
      read(12,1000) zeile
      allcond = " "
      NDUST = 1
      do i=1,imax
        read(12,1000) zeile
        read(zeile,*) dust_nam(NDUST)
        j1 = index(zeile,' ')
        read(zeile(j1+1:),*) trivial(NDUST)
        if (index(zeile,'[l]')>0) then
          j2 = index(zeile,trim(trivial(NDUST)))
     &       + len(trim(trivial(NDUST)))
          read(zeile(j2+1:),*) Tmelt(NDUST)
          trivial(NDUST)=' '
        endif
        read(12,*) dust_rho(NDUST)
        read(12,*) dust_nel(NDUST)
        dmass = 0.d0
        allfound = .true.
        do j=1,dust_nel(NDUST)
          read(12,1030) dust_nu(NDUST,j),name2
          found = .false. 
          do k=1,NELEM
            if (elnam(k).eq.name2) then
              dust_el(NDUST,j) = k
              dmass = dmass + dust_nu(NDUST,j)*mass(k)
              found = .true.
            endif
          enddo
          if (.not.found) then
            print*,trim(dust_nam(NDUST)),name2
            print*,elnam(1:NELEM)
            stop 'Element in dust species not found'
          endif  
          found = .false.
          do k=1,NELM
            if (catm(k).eq.name2) then
              found = .true.
              exit
            endif
          enddo
          if (.not.found) allfound=.false.
        enddo
        found = .false.
        do 
          lastzeile = zeile 
          read(12,1000) zeile
          if (trim(zeile)=='') exit
          if (zeile(1:1)=='#') cycle
          read(zeile,*) fit(NDUST),cfit(NDUST,0:4)
          prec(NDUST) = 0.0
          j1 = index(lastzeile,'+/-')
          j2 = index(lastzeile,':')
          if (j1>0) then
            tmp = lastzeile(j1+3:)
            if (j2>j1) tmp=lastzeile(j1+3:j2-1)            
            read(tmp,*) prec(NDUST)
          endif  
          !print*,trim(tmp),prec(NDUST)
          found = .true.
        enddo
        if (.not.found) then
          print*,"*** syntax error in DustChem.dat, condensate=",
     &         dust_nam(NDUST)
          stop
        endif  
        j1 = index(allcond," "//trim(dust_nam(NDUST)))
        if (j1>0) then
          print*,"*** double condensate in DustChem.dat"
          print*,dust_nam(NDUST)
          stop
        endif  
        if (allfound) then
          dust_mass(NDUST) = dmass
          dust_vol(NDUST) = dmass/dust_rho(NDUST)
          write(*,1060) NDUST,dust_nam(NDUST),dust_rho(NDUST),
     &                  dust_vol(NDUST), (dust_nu(NDUST,j),
     &                  elnam(dust_el(NDUST,j)),j=1,dust_nel(NDUST))
          allcond = " "//trim(allcond)//" "//trim(dust_nam(NDUST))
          NDUST = NDUST+1
        endif
      enddo
      NDUST=NDUST-1
      write(*,*) NDUST," condensed species"
      write(*,*)
      write(*,*) '--- involved elements ---'
      NEPS=0
      elcode(:)=0
      do i=1,NDUST
        do j=1,dust_nel(i)
          name2 = elnam(dust_el(i,j)) 
          do k=1,NELEM
            if (elnam(k).eq.name2) then
              el = k
              exit
            endif
          enddo
          found = .false.           
          do k=1,NEPS
            if (el==elnr(k)) found=.true.
          enddo
          if (.not.found) then
            NEPS = NEPS+1 
            elnr(NEPS) = el
            elcode(el) = NEPS
            write(*,*) elcode(elnr(NEPS)),' ',name2,el
          endif
        enddo
      enddo

      write(*,*) '--- nucleation species ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) NNUC
      allocate(nuc(NNUC),nuc_nam(NNUC))
      do i=1,NNUC
        read(12,*) nuc_nam(i)
        found = .false.
        do j=1,NDUST
          if (trim(nuc_nam(i))//'[s]'==trim(dust_nam(j))) then
            found = .true.
            nuc(i) = j
          endif  
        enddo
        if (found) then
          print'(8x,I4,2x,A6,"->",4x,A10)',
     >          i,nuc_nam(i),dust_nam(nuc(i))
        else
          print*,"*** nucleation species must exist as dust species: ",
     >          trim(nuc_nam(i))
          stop
        endif   
      enddo

      write(*,*) '--- involved gas species ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) NSPECIES
      allocate(spnam(NSPECIES),spnr(NSPECIES),keysp(NSPECIES),
     &         spelem(NSPECIES,10),spstoich(NSPECIES,10),
     &         spzahl(NSPECIES),spmass(NSPECIES))
      do i=1,NSPECIES
        read(12,1020) is_atom,keysp(i),mm,name20
        spnam(i)  = name20
        spnr(i)   = 0
        if (is_atom) then
          spzahl(i) = 1
          name2 = name20(1:2)
          found = .false.
          do j=1,NELEM
            if (name2.eq.elnam(j)) then
              spelem(i,1) = j
              spstoich(i,1) = 1
              spnr(i) = j + 1000
              spmass(i) = mass(j)
              found = .true.
            endif
          enddo
          if (.not.found) stop 'Element of atomic species not found.'
        else
          spnr(i) = stindex(cmol,NMOLE,name20)
          spmass(i) = 0.d0
          spzahl(i) = mm
          do j=1,mm
            read(12,1030) nn,name2
            found = .false. 
            do kk=1,NELEM
              if (elnam(kk).eq.name2) then
                spelem(i,j)   = kk
                spstoich(i,j) = nn 
                spmass(i) = spmass(i) + nn*mass(kk)
                found = .true.
              endif   
            enddo
            if (.not.found) stop 'Element in molecule not found.'
          enddo
        endif
        if (spnr(i).eq.0) stop 'Molecular species not found.'
        write(*,1050) spnam(i),spnr(i),spmass(i)/amu
      enddo

      write(*,*) '--- growth reactions ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) NREAC
      allocate(neduct(NREAC),nprod(NREAC),stick(NREAC),
     >         reac_sp(NREAC,5,2),reac_nu(NREAC,5,2))
      do i=1,NREAC
        do j=1,NELEM
          elem_cons(j) = 0
        enddo  
        read(12,*) neduct(i),nprod(i),stick(i)
        do j=1,neduct(i)
          read(12,1040) reac_nu(i,j,1),name20
          found = .false.
          do kk=1,NSPECIES
            if (name20.eq.spnam(kk)) then
              reac_sp(i,j,1) = kk
              found = .true.
            endif   
          enddo
          if (.not.found) then
            write(*,*) name20 
            stop 'educt in reaction not found'
          endif  
          sp = reac_sp(i,j,1)
          do k=1,spzahl(sp)
            el = spelem(sp,k)  
            elem_cons(el) = elem_cons(el) 
     &                    + reac_nu(i,j,1)*spstoich(sp,k)
          enddo  
        enddo
        read(12,1040) reac_nu(i,1,2),name20
        found = .false.
        do kk=1,NDUST
          if (name20.eq.dust_nam(kk)) then 
            reac_sp(i,1,2) = kk
            sp = reac_sp(i,1,2)
            found = .true.
            do k=1,dust_nel(sp)
              el = dust_el(sp,k)
              elem_cons(el) = elem_cons(el) 
     &                      - reac_nu(i,1,2)*dust_nu(sp,k)
            enddo  
          endif  
        enddo
        if (.not.found) then
          write(*,*) name20
          stop 'dust species in reaction not found'
        endif   
        do j=2,nprod(i)
          read(12,1040) reac_nu(i,j,2),name20
          found = .false.
          do kk=1,NSPECIES
            if (name20.eq.spnam(kk)) then
              reac_sp(i,j,2) = kk
              sp = reac_sp(i,j,2)
              found = .true.
              do k=1,spzahl(sp)
                el = spelem(sp,k)  
                elem_cons(el) = elem_cons(el) 
     &                        - reac_nu(i,j,2)*spstoich(sp,k)
              enddo  
            endif  
          enddo
          if (.not.found) then
            write(*,*) name20
            stop 'product in reaction not found'
          endif   
        enddo
        if (neduct(i).eq.1) then
          write(*,2011) i,stick(i),
     &    (reac_nu(i,j,1),spnam(reac_sp(i,j,1)),j=1,neduct(i)),
     &    reac_nu(i,1,2), dust_nam(reac_sp(i,1,2)), 
     &    (reac_nu(i,j,2),spnam(reac_sp(i,j,2)),j=2,nprod(i))
        elseif (neduct(i).eq.2) then
          write(*,2021) i,stick(i),
     &    (reac_nu(i,j,1),spnam(reac_sp(i,j,1)),j=1,neduct(i)),
     &    reac_nu(i,1,2), dust_nam(reac_sp(i,1,2)), 
     &    (reac_nu(i,j,2),spnam(reac_sp(i,j,2)),j=2,nprod(i))
        elseif (neduct(i).eq.3) then
          write(*,2031) i,stick(i),
     &    (reac_nu(i,j,1),spnam(reac_sp(i,j,1)),j=1,neduct(i)),
     &    reac_nu(i,1,2), dust_nam(reac_sp(i,1,2)), 
     &    (reac_nu(i,j,2),spnam(reac_sp(i,j,2)),j=2,nprod(i))
        else
          write(*,*) i,"neduct>3" 
        endif
        do j=1,NELEM
          if (elem_cons(j).ne.0) then
            write(*,*) elnam(j),elem_cons(j) 
            stop 'element conservation in reaction violated.'
          endif  
        enddo    
      enddo
      close(12)

      Tcorr(:) = -1.d0
      call CHECK_MELTING
      write(*,*)

      !open(unit=1,file='condensates.tex')
      !do i=1,NDUST
      !  limit = ' '
      !  j = index(dust_nam(i),"[l]")
      !  if (Tcorr(i)>0.and.j>0) then
      !    write(limit,'("$>$",I4,"K")') int(Tcorr(i))
      !  else if (Tcorr(i)>0) then
      !    write(limit,'("$<$",I4,"K")') int(Tcorr(i))
      !  endif  
      !  if (prec(i)>0.0) then
      !    write(1,3000)
     &!      i,dust_nam(i),trivial(i),dust_rho(i),
     &!      fit(i),limit,cfit(i,0:4),prec(i)
      !  else  
      !    write(1,3001)
     &!      i,dust_nam(i),trivial(i),dust_rho(i),
     &!      fit(i),limit,cfit(i,0:4)
      !  endif  
      !enddo  
      !close(1)

      RETURN 
 1000 format(a200)
 1010 format(a2)
 1020 format(2(l1,1x),i1,1x,a10)
 1030 format(i2,1x,a2)
 1040 format(i2,1x,a10)
 1050 format(1x,a10,i4,' mass=',0pf7.3," amu")
 1060 format(I4,1x,a20," rhod=",0pf6.3," V0=",1pe11.3,2x,
     &       99(i2,"x",A2,1x))
 2011 format(I3,0pF5.2,1x,1(I2,1x,a8),22x,'->',I2,1x,a10,99(I2,1x,a8))
 2021 format(I3,0pF5.2,1x,2(I2,1x,a8),11x,'->',I2,1x,a10,99(I2,1x,a8))
 2031 format(I3,0pF5.2,1x,3(I2,1x,a8)    ,'->',I2,1x,a10,99(I2,1x,a8))
 3000 format(I3," & ",A20," & ",A25," & ",0pF5.2," & ",I2," & ",A8,
     &       " & ",5(1pE12.5," & "),"$\pm$",0pF4.2,"\\")
 3001 format(I3," & ",A20," & ",A25," & ",0pF5.2," & ",I2," & ",A8,
     &       " & ",5(1pE12.5," & "),9x,"\\")
      end 

***********************************************************************
      SUBROUTINE CHECK_MELTING
***********************************************************************
      use CHEMISTRY,ONLY: NMOLE
      use ELEMENTS,ONLY: NELEM
      use DUST_DATA,ONLY: qp,NDUST,dust_nam,Tmelt,Tcorr,is_liquid
      implicit none
      real*8 :: T
      real(kind=qp) :: nat(NELEM),nmol(NMOLE),Sat(NDUST)
      real(kind=qp) :: old,new,S(NDUST,10000)
      integer :: i,j,k,iT,Ncheck,il,is
      integer :: iliq(NDUST),isol(NDUST)
      character(len=15) :: search

      !--------------------------------------
      ! ***  identify solid/liquid pairs  ***
      !--------------------------------------
      is_liquid(:) = .false.
      Ncheck = 0
      do i=1,NDUST
        k = index(dust_nam(i),'[l]')
        if (k>0) then
          is_liquid(i) = .true. 
          Ncheck = Ncheck+1 
          iliq(Ncheck) = i
          isol(Ncheck) = 0
          search = dust_nam(i)
          search = search(1:k-1)//'[s]'
          do j=1,NDUST
            if (search==dust_nam(j)) then
              isol(Ncheck) = j
            endif
          enddo
          if (isol(Ncheck)==0) then
            print*,"*** liquid without solid "//trim(dust_nam(i))
            stop
          endif  
        endif
      enddo  
      if (Ncheck==0) return

      !-------------------------------
      ! ***  check melting points  ***
      !-------------------------------
      print*
      print*,'auto-correction for spurious liquid <-> solid '//
     &       'phase transitions ...'
      nat = 1.Q-100
      nmol = 1.Q-100
      do iT=100,10000
        T = DBLE(iT) 
        call SUPERSAT(T,nat,nmol,Sat)
        S(:,iT) = Sat(:)
      enddo  
      do i=1,Ncheck
        il = iliq(i)
        is = isol(i)
        do iT=101,10000
          T = DBLE(iT) 
          old = S(is,iT-1)/S(il,iT-1)
          new = S(is,iT)/S(il,iT)
          if (old>1.Q0.and.new<1.Q0) then
            !print'(A15,"-> ",A15,":",2(0pF8.1))',
     &      !     dust_nam(is),dust_nam(il),T,Tmelt(il)
          else if (old<1.Q0.and.new>1.Q0) then
            !print'(A15,"<- ",A15,":",0pF8.1,
     &      !     " false intersection point")',
     &      !     dust_nam(is),dust_nam(il),T
            if (T<Tmelt(il)) then
              Tcorr(il) = 0.5*(T+Tmelt(il))  
              print'(" ... correct ",A15," T <",0pF7.1)',
     &             dust_nam(il),Tcorr(il) 
            else  
              Tcorr(is) = 0.5*(T+Tmelt(il))  !correct solid
              print'(" ... correct ",A15," T >",0pF7.1)',
     &             dust_nam(is),Tcorr(is) 
            endif  
          endif  
        enddo   
      enddo
      do iT=100,10000
        T = DBLE(iT) 
        call SUPERSAT(T,nat,nmol,Sat)
        S(:,iT) = Sat(:)
      enddo  
      print'(26x,"melting point[K]  should be")'
      do i=1,Ncheck
        il = iliq(i)
        is = isol(i)
        do iT=101,10000
          T = DBLE(iT) 
          old = S(is,iT-1)/S(il,iT-1)
          new = S(is,iT)/S(il,iT)
          if (old>1.Q0.and.new<1.Q0) then
            print'(A15,"-> ",A15,":",2(0pF8.1))',
     &           dust_nam(is),dust_nam(il),T,Tmelt(il)
          else if (old<1.Q0.and.new>1.Q0) then
            print'(A15,"<- ",A15,":",0pF8.1,
     &           " false intersection point")',
     &           dust_nam(is),dust_nam(il),T
            stop
          endif  
        enddo   
      enddo
      end
