**********************************************************************
      SUBROUTINE READ_STRUCTURE
**********************************************************************
      use NATURE,ONLY: bar,bk,amu,km,pi
      use PARAMETERS,ONLY: struc_file,Teff,logg,Hp,vzconst,pconst,beta
      use READMODEL,ONLY: Rlay,Tlay,play,rholay,glay,vconvlay,
     >                    zlay,mulay,Difflay,Nlayers
      implicit none
      integer :: elementCount,i
      real :: mixLength,Hplay,pconv,grad,ngas,lmean,Dmicro,Kn,vth,vz
      integer,dimension(1000) :: flag_conv,Z
      logical :: conv

      write(*,*) 
      write(*,*) "reading PHOENIX structure"//trim(struc_file)//" ..."
      write(*,*) "============================================"

      open(42,file=struc_file,status='old')
      do i=1,5 
         read(42,*) 
      enddo
      read(42,100) Teff, logg, mixLength   
      write(*,1000) Teff, logg, mixLength
      read(42,*)
      read(42,*)
      read(42,200) Nlayers
      write(*,*) "Nlayers", Nlayers
      read(42,*)
      read(42,*)
      read(42,200) elementCount
      read(42,*)
      read(42,*)
      read(42,210) (Z(i), i=1,elementCount)
c     write(*,*) (Z(i), i=1,  elementCount)
      read(42,*)
      read(42,*)
      read(42,300) (Rlay(i),i=1,Nlayers)
c     write(*,*) "nach Rlay"
      read(42,*)
      read(42,*)
      read(42,300) (Tlay(i),i=1,Nlayers)
c     write(*,*) "nach Tlay"
      read(42,*)
      read(42,*)
      read(42,300) (play(i),i=1,Nlayers)
c     write(*,*) "nach play"
      read(42,*)
      read(42,*)
      read(42,300) (rholay(i),i=1,Nlayers)
c     write(*,*) "nach rholay"
      read(42,*)
      read(42,*)
      read(42,300) (glay(i),i=1,Nlayers)
c     write(*,*) "nach glay"
      read(42,*)
      read(42,*)
      read(42,300) (vconvlay(i),i=1,Nlayers)
c     write(*,*) "nach vconvlay"
      read(42,*)
      read(42,*)
      read(42,'(8(l2, 1x))') (flag_conv(i),i=1,Nlayers)
      read(42,*)
      read(42,*)
      read(42,300) (mulay(i),i=1,Nlayers)
      read(42,*)
      read(42,*)
c      read(42,300) (epsPhoenix(i), i=1,elementCount)
c      epsH = epsPhoenix(1)
c      do i=1,elementCount
c        if ((Z(i).ge.1).and.(Z(i).le.NELEM)) then
c          eps0(Z(i)) = epsPhoenix(i)/epsH
c        endif  
c      enddo
      close(42)

*     ------------------------------
*     ***  vertical grid points  ***
*     ------------------------------
      do i=1,Nlayers
        zlay(i) = Rlay(i)-Rlay(Nlayers)
      enddo  

*     -----------------------------------------
*     ***  calculate Diffusion coefficient  ***
*     -----------------------------------------
      if (.true.) then
        conv = .true.
        do i=Nlayers,1,-1
          Hplay = bk*Tlay(i)/(glay(i)*mulay(i)*amu)
          if (Tlay(i)>Teff) Hp=Hplay
          if (conv.and.vconvlay(i).gt.0.d0) then
            Difflay(i) = vconvlay(i)*(mixLength*Hplay)  
            !write(*,1200) i,Hplay/km,vconvlay(i)/100.0
          else 
            if (conv) then
              pconv = play(i)
              write(*,'(" convection zone ends at p[bar]=",0pF6.2,
     &                  " (T=",0pF8.2,"K)")') play(i)/bar,Tlay(i)
            endif  
            Difflay(i) = 1.d-99
            conv = .false.
          endif    
          !write(*,*) play(i)/bar,conv,Difflay(i)
        enddo
*       !--- limit log(Diff)-gradient to beta ---
        write(*,*) "Hp[km]=",Hp/km
        do i=Nlayers-1,1,-1
          grad = DLOG(Difflay(i)/Difflay(i+1))/DLOG(play(i)/play(i+1))
          if (grad.gt.beta) then
            Difflay(i) = Difflay(i+1)*DEXP(beta*DLOG(play(i)/play(i+1)))
          endif    
          grad = DLOG(Difflay(i)/Difflay(i+1))/DLOG(play(i)/play(i+1))
          !write(*,1100) i,zlay(i)/km,play(i)/bar,Tlay(i),Difflay(i),grad
        enddo
      else
        do i=1,Nlayers
          vz    = vzconst * MIN(1.d0,play(i)/pconst)**beta
          Hplay = bk*Tlay(i)/(glay(i)*mulay(i)*amu) 
          Difflay(i) = vz*Hplay
        enddo
      endif 

*     -----------------------------
*     ***  add micro diffusion  ***     
*     -----------------------------
      !write(*,'(99(A10))') 'p[bar]','n[cm-3]','vth[km/s]',
     &!                     'l[cm]','Hp[cm]','Dmix','Dmicro'
      do i=Nlayers-1,1,-1
        Hplay  = bk*Tlay(i)/(glay(i)*mulay(i)*amu)
        if (Tlay(i)>Teff) Hp=Hplay
        ngas   = play(i)/bk/Tlay(i)
        lmean  = 1.d0/(2.1E-15*ngas)
        Kn     = lmean/Hp
        vth    = SQRT(8.d0*bk*Tlay(i)/(pi*mulay(i)*amu))
        Dmicro = 1.d0/3.d0*vth*lmean
        !write(*,'(99(1pE10.2))') play(i)/bar,ngas,vth/km,lmean,Hplay,
     >  !        Difflay(i),Dmicro
        Difflay(i) = Difflay(i) + Dmicro
      enddo  

      RETURN
  100 format(f12.3, 1x, f12.3, 1x, f12.3)
  200 format(i5)
  210 format(8(i5, 1x))
  300 format(8(e15.8, 1x)) 
 1000 format(' Teff=',0pF8.3,' logg=',0pF5.2,' mixLengthPara=',0pF5.2)
 1100 format(i4,0pF11.3,1pE11.3,1x,0pF11.2,1x,2(1pe11.3,1x)) 
 1110 format(a4,5(a11,1x)) 
 1200 format(i4,"   Hp[km]=",F8.3,"   vconv[m/s]=",F8.3)
      end
