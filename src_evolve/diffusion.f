************************************************************************
      subroutine DIFFUSION(time,deltat,verbose)
************************************************************************
      use PARAMETERS,ONLY: implicit
      use GRID,ONLY: xlower,xupper,dt_diff_im
      use STRUCT,ONLY: crust_gaseps      
      use ELEMENTS,ONLY: NELEM
      implicit none
      real*8,intent(in) :: time,deltat
      integer,intent(in) :: verbose

      xlower = crust_gaseps

      if (verbose>1) then
        print*
        print*,"entering DIFFUSION ..."
        print*,"======================"
      endif  
      if (implicit.and.deltat>dt_diff_im) then
        call DIFFUSION_IMPLICIT(time,deltat,verbose)
      else
        call DIFFUSION_EXPLICIT(time,deltat,verbose)
      endif
      end


************************************************************************
      subroutine DIFFUSION_EXPLICIT(time0,deltat,verbose)
************************************************************************
      use NATURE,ONLY: pi
      use GRID,ONLY: N=>Npoints,zz,d1l,d1m,d1r,d2l,d2m,d2r,dt_diff_ex,
     >               xlower,xupper
      use PARAMETERS,ONLY: dust_diffuse,bc_low,bc_high,
     >                     influx,outflux,inrate,outrate,vin,vout
      use STRUCT,ONLY: Diff,nHtot,nHeps,crust_Neps
      use DUST_DATA,ONLY: NDUST
      use ELEMENTS,ONLY: elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      implicit none
      real*8,intent(IN) :: time0,deltat
      integer,intent(in) :: verbose
      integer :: i,it,e,el
      real*8,dimension(N) :: xx,rate
      real*8 :: D,nD,d1,d2,d1nD,time,dt
      character :: CR = CHAR(13)
      logical :: IS_NAN

      time = 0.d0
      do it=1,9999999

        dt = MIN(dt_diff_ex,deltat-time)

        do e=1,NELM
          if (e==iel) cycle
          el = elnum(e)
          xx(:) = nHeps(el,:)/nHtot(:) 
          xx(1) = xlower(el)            ! (not sure ...)

          !-------------------------------------------
          ! ***  d/dt(nH*x) = d/dz(nH*Diff*dx/dz)  ***
          !-------------------------------------------
          rate(:) = 0.0
          do i=2,N-1
            nD   = nHtot(i)*Diff(i)  
            d1   = d1l(i)*xx(i-1) + d1m(i)*xx(i) + d1r(i)*xx(i+1)
            d2   = d2l(i)*xx(i-1) + d2m(i)*xx(i) + d2r(i)*xx(i+1)
            d1nD = d1l(i)*nHtot(i-1)*Diff(i-1) 
     >           + d1m(i)*nHtot(i)  *Diff(i) 
     >           + d1r(i)*nHtot(i+1)*Diff(i+1) 
            rate(i) = nD*d2 + d1nD*d1
          enddo

          !----------------------------
          ! ***  explicit timestep  ***
          !----------------------------
          do i=2,N-1
            xx(i) = xx(i) + rate(i)*dt/nHtot(i)
          enddo  

          !------------------------------
          ! ***  boundary conditions  ***
          !------------------------------
          nD = nHtot(1)*Diff(1)  
          if (crust_Neps(el)>0.Q0) then   
            xx(1) = xlower(el)              ! constant concentration
            influx = -nD
     >             *(d1l(1)*xx(1) + d1m(1)*xx(2) + d1r(1)*xx(3))
          else                              
            influx = 0.d0                   ! zero flux 
            xx(1) = (-influx/nD - d1m(1)*xx(2) - d1r(1)*xx(3))
     >             /d1l(1) 
          endif  
          nD = nHtot(N)*Diff(N)  
          if (bc_high==1) then
            xx(N) = xupper(el)              ! const concentration bound.cond. 
            outflux = -nD
     >              *(d1l(N)*xx(N-2) + d1m(N)*xx(N-1) + d1r(N)*xx(N)) 
          else if (bc_high==2) then   
            xx(N) = (-outflux/nD - d1l(N)*xx(N-2) - d1m(N)*xx(N-1))
     >              /d1r(N)
          else if (bc_high==3) then   
            outflux = nHtot(N)*xx(N)*outrate*vout  
            xx(N) = -(d1l(N)*xx(N-2) + d1m(N)*xx(N-1))
     >              /(d1r(N) + outrate*vout/Diff(N))
          endif   

          !----------------------------------------
          ! ***  update crust column densities  ***
          !----------------------------------------
          if (crust_Neps(el)>0.Q0) then   
            print*,elnam(el),REAL(crust_Neps(el)),influx*dt
            crust_Neps(el) = crust_Neps(el) - influx*dt
            if (crust_Neps(el)<0.Q0) then
              print*,elnam(el),"negative crust column density"
              stop
            endif  
          endif  

          !------------------------------------------
          ! ***  map solution on atmosphere grid  ***
          !------------------------------------------
          do i=1,N
            nHeps(el,i) = nHtot(i)*xx(i)
            if (xx(i)<0.d0) then
              !print*,i
              !print*,xx
              !stop "should not occur." 
              if (i==N) then 
                nHeps(el,i) = nHtot(i)*xx(N-1)
              else  
                nHeps(el,i) = nHtot(i)*1.E-50
              endif
            endif    
          enddo 

        enddo  

        time = time + dt
        if (verbose>0) write(*,'(TL10," DIFFUSION:",I8,A,$)') it,CR
        if (time.ge.deltat) exit

      enddo  
      if (verbose>0) print'(" DIFFUSION:",I8,"  time=",1pE11.4,
     >               "  Dt=",1pE11.4)',it,time0,time

      end


************************************************************************
      subroutine DIFFUSION_IMPLICIT(time0,deltat,verbose)
************************************************************************
      use GRID,ONLY: N=>Npoints,zz,d1l,d1m,d1r,d2l,d2m,d2r,
     >               BB,dt_diff_im,xlower,xupper
      use PARAMETERS,ONLY: tfac,bc_low,bc_high,
     >                     influx,outflux,inrate,outrate,vin,vout,
     >                     dust_diffuse
      use STRUCT,ONLY: Diff,nHtot,nHeps,crust_Neps
      use DUST_DATA,ONLY: NDUST
      use ELEMENTS,ONLY: elnam
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      implicit none
      real*8,intent(IN) :: time0,deltat
      integer,intent(in) :: verbose
      integer :: i,j,it,e,el,Nstep
      real*8,dimension(N) :: xx,xnew,rest
      real*8 :: nD,time,dt
      character :: CR = CHAR(13)

      Nstep = deltat/dt_diff_im
      time  = 0.d0
      dt    = deltat/Nstep
      
      do it=1,Nstep

        do e=1,NELM
          if (e==iel) cycle
          el = elnum(e)
          xx(:) = nHeps(el,:)/nHtot(:) 

          !------------------------------
          ! ***  boundary conditions  ***
          !------------------------------
          rest = xx
          nD = nHtot(1)*Diff(1)  
          if (crust_Neps(el)>0.Q0) then   
            xx(1)   = xlower(el)              ! const concentration bound.cond.
            influx  = -nD
     >              *(d1l(1)*xx(1) + d1m(1)*xx(2) + d1r(1)*xx(3))
          else                              
            influx  = 0.d0                    ! zero flux bound.cond.
            rest(1) = -influx/nD/d1l(1) 
          endif
          nD = nHtot(N)*Diff(N)  
          if (bc_high==1) then
            xx(N)   = xupper(el) 
            outflux = -nD
     >              *(d1l(N)*xx(N-2) + d1m(N)*xx(N-1) + d1r(N)*xx(N)) 
          else if (bc_high==2) then
            rest(N) = -outflux/nD/d1r(N)           
          else if (bc_high==3) then
            outflux = outrate*nHtot(N)*xx(N)*vout 
            rest(N) = 0.d0
          endif

          !----------------------------
          ! ***  implicit timestep  ***
          !----------------------------
          do i=1,N
            xnew(i) = 0.d0 
            do j=1,N 
              xnew(i) = xnew(i) + BB(i,j)*rest(j)
            enddo  
          enddo  
          nHeps(el,:) = nHtot(:)*xnew(:)

          !----------------------------------------
          ! ***  update crust column densities  ***
          !----------------------------------------
          if (crust_Neps(el)>0.Q0) then   
            print*,elnam(el),REAL(crust_Neps(el)),influx*dt
            crust_Neps(el) = crust_Neps(el) - influx*dt
            if (crust_Neps(el)<0.Q0) then
              print*,elnam(el),"negative crust column density"
              stop
            endif  
          endif  

        enddo  

        time = time + dt
        if (verbose>0) write(*,'(TL10," DIFFUSION:",I8,A,$)') it,CR

      enddo
      if (verbose>0) print'(" DIFFUSION:",I8,"  time=",1pE11.4,
     >               "  Dt=",1pE11.4)',it,time0,time

      end


