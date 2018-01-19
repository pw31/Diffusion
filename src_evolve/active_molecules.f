************************************************************************
      SUBROUTINE ACTIVE_MOLECULES(nH,Tg,nHeps,verbose)
************************************************************************
      use CHEMISTRY,ONLY: NMOLE,NELEM,NELM,elnum,el,cmol,m_kind,elion
      use ELEMENTS,ONLY: NEPS,eps0,elnr,elnam
      use EXCHANGE,ONLY: inactive,nmol,nel
      use SPECIES,ONLY: NSPECIES,spnr
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in) :: nH,Tg,nHeps(NELEM)
      integer,intent(in) :: verbose
      real(kind=qp) :: eps(NELEM),thres(NELM)
      integer :: i,j,e,sp,Nact
      
      inactive(:) = .false.
      !return
      eps = eps0
      do i=1,NELEM
        eps(e) = nHeps(e)/nH        ! element abundances
      enddo  
      call GGCHEM(nH,Tg,eps,.false.,0)

      inactive(:) = .true.
      do i=1,NELM
        if (i==el) then
          thres(i) = 1.Q-3*nel
        else  
          thres(i) = 1.Q-3*eps(elnum(i))*nH
          !print*,elnam(elnum(i)),real(eps(elnum(i)))
          inactive(elion(i)) = .false.
        endif
      enddo  
      do i=1,NMOLE
        !print*,cmol(i),real(nmol(i)) 
        do j=1,m_kind(0,i)
          e = m_kind(j,i) 
          !print*,elnam(elnum(e)),real(thres(e))
          if (nmol(i)>thres(e)) inactive(i)=.false.
        enddo
      enddo  
      do sp=1,NSPECIES
        i = spnr(sp)
        if (i.gt.1000) cycle    ! atoms
        inactive(i) = .false.
      enddo  

      if (verbose>1) then
        Nact = 0
        do i=1,NMOLE
          if (.not.inactive(i)) Nact=Nact+1 
          !print'(A12,1pE11.3,L2)',cmol(i),nmol(i)/nH,inactive(i) 
        enddo
        print'(I4," active molecules from",I4)',Nact,NMOLE
      endif  

      end
