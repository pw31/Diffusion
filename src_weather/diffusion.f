************************************************************************
      subroutine DIFFUSION(time,deltat,verbose)
************************************************************************
      use PARAMETERS,ONLY: implicit
      implicit none
      real*8,intent(in) :: time,deltat
      integer,intent(in) :: verbose

      if (verbose>1) then
        print*
        print*,"entering DIFFUSION ..."
        print*,"======================"
      endif  
      if (implicit) then
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
      use STRUCT,ONLY: Diff,nHtot,nHeps,rhoLj,rhoL3
      use DUST_DATA,ONLY: NDUST
      use ELEMENTS,ONLY: NEPS
      implicit none
      real*8,intent(IN) :: time0,deltat
      integer,intent(in) :: verbose
      integer :: i,it,el,iwork
      real*8,dimension(N) :: xx,rate
      real*8 :: D,nD,d1,d2,d1nD,time,dt
      character :: CR = CHAR(13)
      logical :: IS_NAN

      iwork = NEPS
      if (dust_diffuse) iwork= NEPS+4+NDUST
      time = 0.d0
      do it=1,9999999

        dt = MIN(dt_diff_ex,deltat-time)

        do el=1,iwork

          if (el<=NEPS) then 
            xx(:) = nHeps(el,:)/nHtot(:) 
          else if (el<=NEPS+4) then
            xx(:) = rhoLj(el-NEPS-1,:)/nHtot(:) 
          else   
            xx(:) = rhoL3(el-NEPS-4,:)/nHtot(:) 
          endif  

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
          if (bc_low==1.or.el>NEPS) then
            xx(1) = xlower(el) 
            influx = -nD
     >             *(d1l(1)*xx(1) + d1m(1)*xx(2) + d1r(1)*xx(3))
          else if (bc_low==2) then   
            xx(1) = (-influx/nD - d1m(1)*xx(2) - d1r(1)*xx(3))
     >             /d1l(1) 
          else if (bc_low==3) then   
            xx(1) = -(d1m(1)*xx(2) + d1r(1)*xx(3))
     >             /(d1l(1) + inrate*vin/Diff(1))
          endif  
          nD = nHtot(N)*Diff(N)  
          if (bc_high==1.or.el>NEPS) then
            xx(N) = xupper(el) 
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

          do i=1,N
            if (el.le.NEPS) then 
              nHeps(el,i) = nHtot(i)*xx(i)
            else if (el.le.NEPS+4) then
              rhoLj(el-NEPS-1,i) = nHtot(i)*xx(i)
            else   
              rhoL3(el-NEPS-4,i) = nHtot(i)*xx(i)
            endif  
            if (xx(i)<0.d0) then
              !print*,i
              !print*,xx
              !stop "should not occur." 
              if (el.le.NEPS) then 
                if (i==N) then 
                  nHeps(el,i) = nHtot(i)*xx(N-1)
                else  
                  nHeps(el,i) = nHtot(i)*1.E-50
                endif
              else if (el.le.NEPS+4) then
                if (i==N) then 
                  rhoLj(el-NEPS-1,i) = nHtot(i)*xx(N-1)
                else  
                  rhoLj(el-NEPS-1,i) = nHtot(i)*1.E-99
                endif
              else   
                if (i==N) then 
                  rhoL3(el-NEPS-4,i) = nHtot(i)*xx(N-1)
                else  
                  rhoL3(el-NEPS-4,i) = nHtot(i)*1.E-99
                endif  
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
     >               BB,dt_diff_im
      use PARAMETERS,ONLY: tfac,bc_low,bc_high,
     >                     influx,outflux,inrate,outrate,vin,vout,
     >                     dust_diffuse
      use STRUCT,ONLY: Diff,nHtot,nHeps,rhoLj,rhoL3
      use DUST_DATA,ONLY: NDUST
      use ELEMENTS,ONLY: NEPS
      implicit none
      real*8,intent(IN) :: time0,deltat
      integer,intent(in) :: verbose
      integer :: i,j,it,el,Nstep
      real*8,dimension(N) :: xx,xnew,rest
      real*8 :: nD,time,dt
      character :: CR = CHAR(13)

      Nstep = deltat/dt_diff_im
      time  = 0.d0
      dt    = deltat/Nstep
      
      do it=1,Nstep

        do el=1,NEPS
          xx(:) = nHeps(el,:)/nHtot(:) 
          !------------------------------
          ! ***  boundary conditions  ***
          !------------------------------
          rest = xx
          nD = nHtot(1)*Diff(1)  
          if (bc_low==1) then
            influx  = -nD
     >              *(d1l(1)*xx(1) + d1m(1)*xx(2) + d1r(1)*xx(3))
          else if (bc_low==2) then
            rest(1) = -influx/nD/d1l(1) 
          else if (bc_low==3) then
            influx  = inrate*nHtot(1)*xx(1)*vin 
            rest(1) = 0.d0
          endif
          nD = nHtot(N)*Diff(N)  
          if (bc_high==1) then
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
        enddo  

        time = time + dt
        if (verbose>0) write(*,'(TL10," DIFFUSION:",I8,A,$)') it,CR

      enddo
      if (verbose>0) print'(" DIFFUSION:",I8,"  time=",1pE11.4,
     >               "  Dt=",1pE11.4)',it,time0,time

      end


