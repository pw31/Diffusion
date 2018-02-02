*********************************************************************
      SUBROUTINE JSTAR(T,Sat,nat,nmol,Jst,Nst)
*********************************************************************
      use NUCLEATION,ONLY: NNUC,nuc,nuc_nam
      use NATURE,ONLY: pi
      use ELEMENTS,ONLY: NELEM
      use CHEMISTRY,ONLY: NMOLE,cmol
      use DUST_DATA,ONLY: NDUST,qp,dust_vol
      use EXCHANGE,ONLY: C
      implicit none
      real*8,intent(in) :: T
      real(kind=qp),intent(in) :: nat(NELEM),nmol(NMOLE)
      real(kind=qp),intent(in) :: Sat(NDUST)
      real*8,intent(out):: Nst(NNUC),Jst(NNUC)
      real*8 :: nTiO2,nSiO,SS,Jstern,Nstern
      real*8 :: nC1,nC2,nC2H,nC2H2,nC3,nKCl,r1
      integer :: stindex,i
      integer,save :: TiO2,C2,C3,C2H,C2H2,SiO,KCl
      logical,save :: firstCall=.true.
!$omp threadprivate(firstCall,TiO2,C2,C3,C2H,C2H2,SiO,KCl)

      if (firstCall) then
        TiO2 = STINDEX(CMOL,NMOLE,'TIO2     ')
        SiO  = STINDEX(CMOL,NMOLE,'SIO      ')
        C2   = STINDEX(CMOL,NMOLE,'C2       ') 
        C3   = STINDEX(CMOL,NMOLE,'C3       ')
        C2H  = STINDEX(CMOL,NMOLE,'C2H      ')
        C2H2 = STINDEX(CMOL,NMOLE,'C2H2     ')
	KCl  = STINDEX(CMOL,NMOLE,'KCL      ')
        firstCall=.false.
      endif    

      do i=1,NNUC
        if (nuc_nam(i).eq.'TiO2') then
          nTiO2 = nmol(TiO2)
          SS    = Sat(nuc(i))
          CALL CLASS_NUC_TIO2(T,SS,nTiO2,Jstern,Nstern)
          !write(*,*) "TiO2 nuc:",nTiO2,SS,Jstern,J_is_zero
        else if (nuc_nam(i).eq.'SiO') then
          nSiO  = nmol(SiO)
          SS    = Sat(nuc(i))
          CALL CLASS_NUC_SIO(T,SS,nSIO,Jstern,Nstern)
          Jst(i) = Jstern
          !write(*,*) "SiO nuc in NUCLEATION.f"
        else if (nuc_nam(i).eq.'C') then
          nC1   = nat(C)
          nC2   = nmol(C2)
          nC3   = nmol(C3)
          nC2H  = nmol(C2H)
          nC2H2 = nmol(C2H2)
          SS    = Sat(nuc(i))
          CALL CLASS_NUC_C(T,SS,nC1,nC2,nC2H,nC2H2,nC3,Jstern,Nstern)
          !write(*,*) "C nuc:",nC1,SS,Jstern,J_is_zero
	else if (nuc_nam(i).eq.'KCl') then
          nKCl = nmol(KCl)
          SS = Sat(nuc(i))
          CALL CLASS_NUC_KCL(T,SS,nKCl,Jstern,Nstern)
        else
          write(*,*) 'unknown nucleation species',nuc_nam(i)
        endif
        Jst(i) = Jstern
        Nst(i) = Nstern
        !r1 = (dust_vol(nuc(i))/(4.d0*pi)*3.d0)**(1.d0/3.d0)
        !print*,nuc_nam(i),4*pi*r1**2
      enddo

      RETURN 
      END 
