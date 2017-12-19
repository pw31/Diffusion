************************************************************************
      SUBROUTINE SETTLING(time,deltat,verbose)
************************************************************************
      use NATURE,ONLY: pi,bk,amu,mic
      use PARAMETERS,ONLY: logg
      use ELEMENTS,ONLY: NEPS,muH
      use DUST_DATA,ONLY: NDUST,dust_rho,dust_vol
      use GRID,ONLY: Npoints,zz
      use STRUCT,ONLY: nHtot,Temp,mu,nHeps,rhoLj,rhoL3
      implicit none
      integer,intent(in) :: verbose
      real*8,intent(in) :: time
      real*8,intent(inout) :: deltat
      real*8 :: g,xi,rho,Tg,rhod,cT,bsum,vdrift,delz,dt,tpass,epsd
      real*8 :: LL(0:4),bmix(NDUST),jj(4+NDUST,Npoints)
      real*8 :: jin(4+NDUST),jout(4+NDUST),CLOSURE
      integer :: ip,i,imin

      if (verbose>1) then
        print*
        print*,"entering SETTLING ..."
        print*,"====================="
      endif  

      g  = 10.d0**logg
      xi = DSQRT(pi)/2.d0 * (3.d0/(4.d0*pi))**(1.d0/3.d0) * g
      dt = 1.d+99
      imin = 0

      jj(:,:) = 0.d0
      do ip=1,Npoints
        rho = nHtot(ip)*muH
        Tg  = Temp(ip)
        cT  = DSQRT(2.d0*bk*Tg/mu(ip))
        epsd = rhoLj(3,ip)/dust_Vol(1)/nHtot(ip)  ! cm^3/cm^3 / cm3 / cm-3   
        if (epsd<1.E-15) then                     ! dust-free case
          if (ip>1) then
            rhod = 3.0
            vdrift = xi*rhod/cT/rho*(1.E-2*mic)
            delz = ABS(zz(ip-1)-zz(ip)) 
            dt = MIN(dt,0.6*delz/vdrift)
            !print'(I3,99(1pE12.3))',ip,Tg,epsd,vdrift,dt
          endif  
          cycle
        endif  

        !--- closure condition ---
        LL(0:3) = rhoLj(0:3,ip)/rho               ! [cm^j/g]
        LL(4)   = CLOSURE(ip,LL(0),LL(1),LL(2),LL(3),verbose)
        rhoLj(1,ip) = rho*LL(1)
        rhoLj(2,ip) = rho*LL(2)

        !--- dust material density ---
        bmix(1:NDUST) = rhoL3(1:NDUST,ip)    
        bsum = 0.d0                         
        rhod = 0.d0                               
        do i=1,NDUST
          bsum = bsum + bmix(i)
          rhod = rhod + bmix(i)*dust_rho(i)
        enddo
        if (bsum>0.d0) then
          bmix = bmix/bsum
          rhod = rhod/bsum
        else
          stop "*** should not happen."
        endif  
        do i=0,3
          jj(i+1,ip) = xi*rhod/cT*LL(i+1) 
          vdrift     = jj(i+1,ip)/rhoLj(i,ip)
          !vdrift     = xi*rhod/cT/rho*(LL(i+1)/LL(i)) 
          !jj(i+1,ip) = vdrift*rhoLj(i,ip)
          if (ip>1) then
            delz = ABS(zz(ip-1)-zz(ip)) 
            tpass = 0.6*delz/vdrift
            !print'(I4,I4,99(1pE10.2))',ip,i,delz,xi,rhod,cT,vdrift
            if (tpass<dt) then
              dt = tpass                          ! limit timestep
              imin = ip
              !print'(2(I3),99(1pE12.3))',ip,i,epsd,LL,vdrift
            endif  
          endif  
        enddo
        vdrift = xi*rhod/cT*LL(4)/rhoLj(3,ip) 
        !vdrift = xi*rhod/cT/rho*(LL(4)/LL(3))
        do i=1,NDUST
          jj(4+i,ip) = vdrift*rhoL3(i,ip)
        enddo
      enddo  
      deltat = min(2.0*deltat,dt)
      print'(" SETTLING: ",I3,"  Dt =",2(1pE10.3))',imin,dt,deltat

      do ip=1,Npoints
        jout(:) = jj(:,ip)                        ! upwind scheme 
        jin(:)  = 0.d0
        if (ip<Npoints) jin(:)=jj(:,ip+1)         ! upwind scheme 

        if (ip==1) then
          delz = ABS(zz(2)-zz(1))
        else if (ip==Npoints) then 
          delz = ABS(zz(Npoints)-zz(Npoints-1))
        else
          delz = 0.5*ABS(zz(ip+1)-zz(ip-1))
        endif   

        rhoLj(0:3,ip) = rhoLj(0:3,ip) 
     >     + (jin(1:4) - jout(1:4))*deltat/delz
        rhoL3(1:NDUST,ip) = rhoL3(1:NDUST,ip) 
     >     + (jin(5:4+NDUST) - jout(5:4+NDUST))*deltat/delz

      enddo

      end
