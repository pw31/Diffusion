************************************************************************
      subroutine INITIAL_CONDITIONS(nout,time,dt)
************************************************************************
      use NATURE,ONLY: pi
      use PARAMETERS,ONLY: init,model_name
      use GRID,ONLY: N=>Npoints,xlower,xupper
      use STRUCT,ONLY: nHtot,nHeps,crust_depth,
     >                 crust_beta,crust_Ncond,crust_Neps,crust_gaseps
      use DUST_DATA,ONLY: NDUST
      use ELEMENTS,ONLY: NELEM,eps_solar
      implicit none
      integer,intent(out) :: nout
      real*8,intent(out) :: time,dt
      integer :: i,el
      logical :: ex

      allocate(nHeps(NELEM,-2:N))
      allocate(xlower(NELEM),xupper(NELEM))
      inquire(file=trim(model_name)//"/restart.dat",exist=ex)
      if (init==-1.and.(.not.ex)) init=0

      if (init==-1) then

        open(70,file=trim(model_name)//"/restart.dat",
     >       form="unformatted",status="old")
        read(70) nout,time,dt
        read(70) nHeps
        read(70) crust_depth
        read(70) crust_beta
        read(70) crust_Ncond
        read(70) crust_Neps
        read(70) crust_gaseps
        close(70)
        print*
        print*,"restart from file:  num,time,dt =",nout,time,dt

      else  

        do i=1,N
          if (init==2) then  
            nHeps(:,i) = nHtot(i)*eps_solar(:)     ! solar abundances  
          else if (init==1) then
            nHeps(:,i) = nHtot(i)*crust_gaseps(:)  ! gas abundances over crust
          else if (init==0) then
            nHeps(:,i) = nHtot(i)*1.d-50           ! empty        
          else
            print*,"*** init=",init," not recognised."
            stop
          endif  
        enddo
        do i=-2,2    ! guard cells -2,-1 and first three cells 0,1,2
          nHeps(:,i) = nHtot(i)*crust_gaseps(:)
        enddo  
        print*
        print*,"start from init =",init

      endif  

      xlower(:) = crust_gaseps(:)
      xupper(:) = nHeps(:,N)/nHtot(N)

      end
