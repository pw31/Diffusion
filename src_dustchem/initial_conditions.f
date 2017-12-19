************************************************************************
      subroutine INITIAL_CONDITIONS(nout,time,dt)
************************************************************************
      use NATURE,ONLY: pi
      use PARAMETERS,ONLY: init,model_name
      use GRID,ONLY: N=>Npoints,zz,xlower,xupper
      use STRUCT,ONLY: nHtot,Diff,nHeps,rhoLj,rhoL3
      use DUST_DATA,ONLY: NDUST
      use ELEMENTS,ONLY: NEPS,elnr,eps0
      use EXCHANGE,ONLY: C,O,S,Cl
      implicit none
      integer,intent(out) :: nout
      real*8,intent(out) :: time,dt
      real*8 :: esolar(NEPS),eempty(NEPS)
      integer :: i,el
      logical :: ex

      allocate(nHeps(NEPS,N))
      allocate(rhoLj(0:3,N))
      allocate(rhoL3(NDUST,N))
      allocate(xlower(4+NDUST+NEPS),xupper(4+NDUST+NEPS))
      inquire(file=trim(model_name)//"/restart.dat",exist=ex)
      if (init==-1.and.(.not.ex)) init=0

      if (init==-1) then

        open(70,file=trim(model_name)//"/restart.dat",
     >       form="unformatted",status="old")
        read(70) nout,time,dt
        read(70) nHeps
        read(70) rhoLj
        read(70) rhoL3
        close(70)
        print*
        print*,"restart from file:  num,time,dt =",nout,time,dt

      else  

        do i=1,NEPS
          el = elnr(i) 
          esolar(i) = eps0(el)
          eempty(i) = 1.d-50
          if (el==C)  eempty(i)=eps0(C) 
          if (el==O)  eempty(i)=0.7*eps0(O) 
          if (el==S)  eempty(i)=eps0(S) 
          if (el==Cl) eempty(i)=eps0(Cl) 
        enddo  
        rhoLj = 0.d0          ! dust-free
        rhoL3 = 0.d0          ! dust-free
        do i=1,N
          if (init==0) then  
            if (i.eq.1) then 
              nHeps(1:NEPS,i) = nHtot(i)*esolar(1:NEPS)  
            else
              nHeps(1:NEPS,i) = nHtot(i)*eempty(1:NEPS) 
            endif  
          else if (init==1) then
            nHeps(1:NEPS,i) = nHtot(i)*esolar(1:NEPS) 
          else
            print*,"*** init=",init," not recognised."
            stop
          endif  
        enddo
        print*
        print*,"start from init =",init
      endif  

      xlower(1:NEPS) = nHeps(1:NEPS,1)/nHtot(1)
      xupper(1:NEPS) = nHeps(1:NEPS,N)/nHtot(N)
      xlower(NEPS+1:NEPS+4+NDUST) = 0.d0 
      xupper(NEPS+1:NEPS+4+NDUST) = 0.d0 

      end
