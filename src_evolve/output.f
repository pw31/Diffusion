***********************************************************************
      SUBROUTINE OUTPUT(num,time,dt)
***********************************************************************
      use PARAMETERS,ONLY: model_name,logg,verbose
      use NATURE,ONLY: bk,bar,amu,mel,pi,mic
      use CHEMISTRY,ONLY: NELM,NMOLE,elnum,cmol,catm,el,charge,molmass
      use DUST_DATA,ONLY: NDUST,dust_nam,dust_mass,dust_Vol,
     &                    dust_nel,dust_el,dust_nu,dust_rho
      use EXCHANGE,ONLY: ipoint,nel,nat,nion,nmol,Jst,chi
      use GRID,ONLY: Npoints,zz,d1l,d1m,d1r,zweight
      use STRUCT,ONLY: Temp,press,rho,nHtot,Diff,nHeps,crust_gaseps,
     &                 crust_depth,crust_Ncond,crust_Neps,crust_beta
      use NUCLEATION,ONLY: NNUC,nuc,nuc_nam
      use ELEMENTS,ONLY: NELEM,elnr,elcode,elnam,eps0,mass,muH
      implicit none
      integer,intent(in) :: num
      real*8,intent(in) :: time,dt
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST),out(NDUST)
      real*8 :: pp,Tg,nH,nges,kT,mu,sumn,sumnm,dz,Ncol(NELEM)
      integer :: i,j,ip,dk,e,NOUT
      character(len=200) :: line,filename
      character(len=20) :: name,short_name(NDUST)
      character(len=1) :: char1
      character(len=8) :: cout
      logical :: hasW,ex

#ifdef IFORT
      inquire(directory=trim(model_name),exist=ex)
#else
      inquire(file=trim(model_name)//"/.",exist=ex)
#endif
      if (.not.ex) then
        call SYSTEM("mkdir "//trim(model_name)) 
      endif
  
      !---------------------------------
      ! ***  write crust properties  ***
      !---------------------------------
      inquire(file=trim(model_name)//'/history1.out',exist=ex)
      if ((num==0).or.(.not.ex)) then
        open(unit=12,file=trim(model_name)//'/history1.out',
     >       status='replace') 
        write(12,'(a6,99(a14))') '#','time[s]','dt[s]','depth[cm]',
     >            elnam(elnum(1:el-1)),elnam(elnum(el+1:NELM))
      else   
        open(unit=12,file=trim(model_name)//'/history1.out',
     >       position='append')
      endif   
      write(12,'(i6,99(1pE14.6))') num,time,dt,crust_depth,
     >            crust_Neps(elnum(1:el-1)),crust_Neps(elnum(el+1:NELM))
      close(12)
      inquire(file=trim(model_name)//'/history2.out',exist=ex)
      if ((num==0).or.(.not.ex)) then
        open(unit=12,file=trim(model_name)//'/history2.out',
     >       status='replace') 
        write(12,'(a6,999(a20))') '#','time[s]','dt[s]','depth[cm]',
     >            (trim(dust_nam(i)),i=1,NDUST)
      else   
        open(unit=12,file=trim(model_name)//'/history2.out',
     >       position='append')
      endif   
      write(12,'(i6,999(1pE20.6))') num,time,dt,crust_depth,
     >            crust_Ncond(1:NDUST)
      close(12)
      
      !---------------------------------------
      ! ***  write total column densities  ***
      !---------------------------------------
      inquire(file=trim(model_name)//'/check.out',exist=ex)
      if ((num==0).or.(.not.ex)) then
        open(unit=12,file=trim(model_name)//'/check.out',
     >       status='replace') 
        write(12,'(a6,99(a14))') '#','time[s]','dt[s]','depth[cm]',
     >            elnam(elnum(1:el-1)),elnam(elnum(el+1:NELM))
      else   
        open(unit=12,file=trim(model_name)//'/check.out',
     >       position='append')
      endif   
      Ncol(:) = crust_Neps(:)
      do ip=0,Npoints
        Ncol(:) = Ncol(:) + nHeps(:,ip)*zweight(ip)
      enddo  
      write(12,'(i6,99(1pE14.6))') num,time,dt,crust_depth,
     >            Ncol(elnum(1:el-1)),Ncol(elnum(el+1:NELM))
      close(12)
      
      !---------------------------------------------------------------
      ! ***  write crust and atmospheric structure to output file  ***
      !---------------------------------------------------------------
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
      write(70,*) NOUT,NMOLE,NDUST,Npoints+1
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

      do ip=0,Npoints
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
      if (verbose>=2) read'(A1)',char1

 1000 format(4(' eps(',a2,') = ',1pD8.2))
 1010 format(A4,0pF8.2,3(a6,1pE9.2),1(a11,1pE9.2))
 2000 format(9999(1x,A19))
 2010 format(0pF20.6,3(1pE20.6),9999(0pF20.7))
 3000 format(A20,2(1pE12.4))
      end  

