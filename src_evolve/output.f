***********************************************************************
      SUBROUTINE OUTPUT(num,time,dt)
***********************************************************************
      use PARAMETERS,ONLY: model_name,logg
      use NATURE,ONLY: bk,bar,amu,mel,pi,mic
      use CHEMISTRY,ONLY: NELM,NMOLE,elnum,cmol,catm,el,charge,molmass
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_mass,dust_Vol,
     &                    dust_nel,dust_el,dust_nu,dust_rho
      use EXCHANGE,ONLY: ipoint,nel,nat,nion,nmol,Jst,chi
      use GRID,ONLY: Npoints,zz,d1l,d1m,d1r
      use STRUCT,ONLY: Temp,press,rho,nHtot,Diff,nHeps,crust_gaseps,
     &                 crust_depth,crust_Ncond,crust_Neps,crust_beta
      use NUCLEATION,ONLY: NNUC,nuc,nuc_nam
      use ELEMENTS,ONLY: NELEM,elnr,elcode,elnam,eps0,mass,muH
      implicit none
      integer,intent(in) :: num
      real*8,intent(in) :: time,dt
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST),out(NDUST)
      real*8 :: pp,Tg,nH,nges,kT,mu,sumn,sumnm
      integer :: i,j,ip,dk,e,NOUT
      character(len=200) :: line,filename
      character(len=20) :: name,short_name(NDUST)
      character(len=1) :: char
      character(len=8) :: cout
      logical :: hasW,ex
      integer :: verbose=0

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
      open(unit=70,file=trim(model_name)//'/structure_'//cout//'.dat',
     &     status='replace')
      write(70,*) 't=',time
      write(70,*) NOUT,NMOLE,NDUST,Npoints
      write(70,*) crust_depth
      do i=1,NDUST
        write(70,3000) trim(short_name(i)),crust_Ncond(i),crust_beta(i) 
      enddo  
      write(70,2000) 'Tg','nHges','pges','Dmix','el',
     &               (trim(elnam(elnum(j))),j=1,el-1),
     &               (trim(elnam(elnum(j))),j=el+1,NELM),
     &               (trim(cmol(i)),i=1,NMOLE),
     &               ('S'//trim(short_name(i)),i=1,NDUST),
     &               ('eps'//trim(elnam(elnum(j))),j=1,el-1),
     &               ('eps'//trim(elnam(elnum(j))),j=el+1,NELM)

      do ip=1,Npoints
        ipoint = ip

        !--- temperature and density ---
        nH = nHtot(ip)
        Tg = Temp(ip) 
        kT = bk*Tg

        !--- compute chemistry, supersaturation, nucleation ---
        eps = eps0
        do i=1,NELEM
          eps(i) = nHeps(i,ip)/nH
        enddo         
        call GGCHEM(nH,Tg,eps,.false.,0)
        call SUPERSAT(Tg,nat,nmol,Sat)

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
        
        write(70,2010) Tg,nH,pp,Diff(ip),
     &       LOG10(MAX(1.Q-300, nel)),   
     &      (LOG10(MAX(1.Q-300, nat(elnum(j)))),j=1,el-1),
     &      (LOG10(MAX(1.Q-300, nat(elnum(j)))),j=el+1,NELM),
     &      (LOG10(MAX(1.Q-300, nmol(j))),j=1,NMOLE),
     &      (LOG10(MAX(1.Q-300, Sat(j))),j=1,NDUST),
     &      (LOG10(eps(elnum(j))),j=1,el-1),
     &      (LOG10(eps(elnum(j))),j=el+1,NELM)
        
      enddo  
      close(70)

      open(70,file=trim(model_name)//'/restart.dat',
     &     form="unformatted",status="replace")
      write(70) num,time,dt
      write(70) nHeps
      write(70) crust_depth
      write(70) crust_beta
      write(70) crust_Ncond
      write(70) crust_Neps
      write(70) crust_gaseps
      close(70)

 1000 format(4(' eps(',a2,') = ',1pD8.2))
 1010 format(A4,0pF8.2,3(a6,1pE9.2),1(a11,1pE9.2))
 2000 format(9999(1x,A19))
 2010 format(0pF20.6,3(1pE20.6),9999(0pF20.7))
 3000 format(A20,2(1pE12.4))
      end  

