************************************************************************
      subroutine INITIAL_CONDITIONS(nout,time,tnext,dt)
************************************************************************
      use NATURE,ONLY: pi
      use PARAMETERS,ONLY: gas_kind,model_name
      use GRID,ONLY: N=>Npoints,xlower,xupper
      use STRUCT,ONLY: nHtot,nHeps,crust_depth,
     >                 crust_beta,crust_Ncond,crust_Neps,crust_gaseps
      use DUST_DATA,ONLY: NDUST
      use ELEMENTS,ONLY: NELEM,eps_solar,eps_meteor,eps_crust
      use EXCHANGE,ONLY: H,He
      use JEANS_ESCAPE,ONLY: EXTRA
      implicit none
      integer,intent(inout) :: nout
      real*8,intent(inout) :: time,tnext,dt
      integer :: i,el,Nsolve,indep(NELEM)
      logical :: ex

      allocate(nHeps(NELEM+EXTRA,-2:N))
      allocate(xlower(NELEM),xupper(NELEM))
      inquire(file=trim(model_name)//"/restart.dat",exist=ex)
      if (gas_kind==-1.and.(.not.ex)) then
        stop "*** no restart file."
      endif  

      if (gas_kind==-1) then

        open(70,file=trim(model_name)//"/restart.dat",
     >       form="unformatted",status="old")
        read(70) nout,time,tnext,dt
        read(70) nHeps
        read(70) crust_depth
        read(70) crust_beta
        read(70) crust_Ncond
        read(70) crust_Neps
        read(70) crust_gaseps
        close(70)
        print*
        print*,"restart from file:  num,time,dt =",nout,time,dt
        call WARM_UP(time)
        call UPDATE_CRUST(Nsolve,indep)

      else  

        do i=1,N
          if (gas_kind==0) then
            nHeps(:,i) = nHtot(i)*1.d-50           ! empty        
            nHeps(H,i) = nHtot(i)*eps_solar(H)
            nHeps(He,i)= nHtot(i)*eps_solar(He)
          else if (gas_kind==1) then  
            nHeps(:,i) = nHtot(i)*eps_solar(:)     ! solar abundances  
          else if (gas_kind==2) then
            nHeps(:,i) = nHtot(i)*eps_meteor(:)    ! meteoritic abundances
          else if (gas_kind==3) then
            nHeps(:,i) = nHtot(i)*eps_crust(:)     ! Earth crust abundances
          else if (gas_kind==4) then
            nHeps(:,i) = nHtot(i)*crust_gaseps(:)  ! gas abundances over crust
          else
            print*,"*** gas_kind=",gas_kind," not supported."
            stop
          endif  
        enddo
        do i=-2,1    ! guard cells -2,-1 and first three cells 0,1,2
          nHeps(:,i) = nHtot(i)*crust_gaseps(:)
        enddo  
        print*,"INITIAL_CONDITIONS gas_kind =",gas_kind
        print*

      endif  

      xlower(:) = crust_gaseps(:)
      xupper(:) = nHeps(:,N)/nHtot(N)

      end
