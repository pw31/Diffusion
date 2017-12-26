***********************************************************************
      SUBROUTINE OUTPUT(num,time,dt)
***********************************************************************
      use PARAMETERS,ONLY: model_name
      use NATURE,ONLY: bk,bar,amu,mel,pi,mic
      use CHEMISTRY,ONLY: NELM,NMOLE,elnum,cmol,catm,el,charge,molmass
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_mass,dust_Vol
      use EXCHANGE,ONLY: ipoint,nel,nat,nion,nmol,Jst,chi
      use GRID,ONLY: Npoints,zz
      use STRUCT,ONLY: Temp,press,rho,nHtot,Diff,nHeps,rhoLj,rhoL3
      use NUCLEATION,ONLY: NNUC,nuc,nuc_nam
      use ELEMENTS,ONLY: NELEM,NEPS,elnr,elnam,eps0,mass,muH
      implicit none
      integer,intent(in) :: num
      real*8,intent(in) :: time,dt
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST),out(NDUST)
      real*8 :: pp,Tg,nH,nges,kT,mu,sumn,sumnm
      real*8 :: rhog,dustV,rhod,amean,amax
      real*8 :: rhoL(0:3),Nst(NNUC),bmix(NDUST),effSat(NDUST)
      integer :: i,j,ip,e,NOUT
      character(len=200) :: line,filename
      character(len=20) :: name,short_name(NDUST)
      character(len=1) :: char
      character(len=8) :: cout
      logical :: hasW,ex
      integer :: verbose=0

      if (.not.allocated(nmol)) then
        allocate(nmol(NMOLE),Jst(NNUC),chi(NDUST))
      endif
  
#ifdef IFORT
      inquire(directory=trim(model_name),exist=ex)
#else
      inquire(file=trim(model_name)//"/.",exist=ex)
#endif
      if (.not.ex) then
        call SYSTEM("mkdir "//trim(model_name)) 
      endif  

      !----------------------------
      ! ***  open output files  ***
      !----------------------------
      do i=1,NDUST
        name = dust_nam(i) 
        j=index(name,"[s]")
        short_name(i) = name
        if (j>0) short_name(i)=name(1:j-1)
      enddo
      NOUT = NELM
      if (charge) NOUT=NOUT-1
      write(cout,'(I8.8)') num
      open(unit=70,file=trim(model_name)//'/weather_'//cout//'.dat',
     &     status='replace')
      write(70,*) 't=',time
      write(70,*) NOUT,NMOLE,NDUST,NEPS,NNUC,Npoints
      write(70,2000) 'Tg','nHges','pges','Dmix','el',
     &               (trim(elnam(elnum(j))),j=1,el-1),
     &               (trim(elnam(elnum(j))),j=el+1,NELM),
     &               (trim(cmol(i)),i=1,NMOLE),
     &               ('S'//trim(short_name(i)),i=1,NDUST),
     &               ('n'//trim(short_name(i)),i=1,NDUST),
     &               ('eps'//trim(elnam(elnum(j))),j=1,el-1),
     &               ('eps'//trim(elnam(elnum(j))),j=el+1,NELM),
     &               'dust/gas','dustVol/H','<a>[mic]',
     &               ('Jstar('//trim(nuc_nam(j))//')',j=1,NNUC),
     &               ('Nstar('//trim(nuc_nam(j))//')',j=1,NNUC)

      amax = 0.d0
      do ip=1,Npoints
        ipoint = ip

        !--- temperature and density ---
        nH = nHtot(ip)
        Tg = Temp(ip) 
        kT = bk*Tg

        !--- dust densities and volume mixing ratios ---
        rhoL(0:3) = rhoLj(0:3,ip)
        do j=1,NDUST
          eldust(j) = rhoL3(j,ip)/dust_Vol(j)/nH
        enddo 
        if (rhoL(0)>0.d0) then
          !amean  = (3.d0/(4.d0*pi))**(1.d0/3.d0) * rhoL(1)/rhoL(0)
          amean   = (3.d0/(4.d0*pi)*rhoL(3)/rhoL(0))**(1.d0/3.d0)
          bmix(:) = rhoL3(:,ip)/rhoL(3)
        else
          amean   = 0.d0
          bmix(:) = 1.d0
        endif  

        !--- compute chemistry, supersaturation, nucleation ---
        eps = eps0
        do i=1,NEPS
          e = elnr(i) 
          eps(e) = nHeps(i,ip)/nH
        enddo         
        !call GGCHEM(nH,Tg,eps,.false.,0)
        !call SUPERSAT(Tg,nat,nmol,Sat)
        call EFF_SUPERSAT(nH,Tg,eps,bmix,Sat,effSat)
        call JSTAR(Tg,Sat,nat,nmol,Jst,Nst)

        !--- compute pressure and mean molecular mass ---
        sumn  = nel
        sumnm = nel*mel
        do i=1,NELEM
          if (nat(i)==0.Q0) cycle
          sumn  = sumn  + nat(i)
          sumnm = sumnm + nat(i)*mass(i)
          !write(*,1000) trim(elnam(i)),nat(i),mass(i)/amu
        enddo  
        do i=1,NMOLE
          sumn  = sumn  + nmol(i)
          sumnm = sumnm + nmol(i)*molmass(i)
          !write(*,1000) trim(cmol(i)),nmol(i),molmass(i)/amu
        enddo
        mu = sumnm/sumn
        pp = sumn*bk*Tg
        
        !--- compute dust/gas mass ratio ---
        rhog  = nH*muH
        rhod  = 0.0
        dustV = 0.0
        do j=1,NDUST
          rhod  = rhod  + nH*eldust(j)*dust_mass(j)
          dustV = dustV + eldust(j)*dust_Vol(j)
          out(j) = LOG10(MIN(1.Q+300,MAX(1.Q-300,effSat(j))))
        enddo  

        write(70,2010) Tg,nH,pp,Diff(ip),
     &       LOG10(MAX(1.Q-300, nel)),   
     &      (LOG10(MAX(1.Q-300, nat(elnum(j)))),j=1,el-1),
     &      (LOG10(MAX(1.Q-300, nat(elnum(j)))),j=el+1,NELM),
     &      (LOG10(MAX(1.Q-300, nmol(j))),j=1,NMOLE),
     &      (out(j),j=1,NDUST),
     &      (LOG10(MAX(1.Q-300, eldust(j))),j=1,NDUST),
     &      (LOG10(eps(elnum(j))),j=1,el-1),
     &      (LOG10(eps(elnum(j))),j=el+1,NELM),
     &       LOG10(MAX(1.Q-300, rhod/rhog)),
     &       LOG10(MAX(1.Q-300, dustV)),
     &       amean/mic,
     &      (LOG10(MAX(1.Q-300, Jst(j))),j=1,NNUC), 
     &      (MIN(999999.99999,Nst(j)),j=1,NNUC)
        
        amax = max(amax,amean)
      enddo  
      close(70)

      print'("amax[mic] =",1pE13.5)',amax

      open(70,file=trim(model_name)//'/restart.dat',
     &     form="unformatted",status="replace")
      write(70) num,time,dt
      write(70) nHeps
      write(70) rhoLj
      write(70) rhoL3
      close(70)

 1000 format(4(' eps(',a2,') = ',1pD8.2))
 1010 format(A4,0pF8.2,3(a6,1pE9.2),1(a11,1pE9.2))
 2000 format(9999(1x,A19))
 2010 format(0pF20.6,3(1pE20.6),9999(0pF20.7))
      end  

