************************************************************************
      REAL*8 FUNCTION CLOSURE(ip,L0,L1,L2,L3,verbose)
************************************************************************
      use NATURE,ONLY: mic
      implicit none
      integer,intent(in) :: ip,verbose
      real*8,intent(inout) :: L0,L1,L2,L3
      real*8 :: n1,n2,a1,a2
      real*8 :: aa,bb,cc,pp,qq,VIETA1,VIETA2,VIETA3,VIETA4
      integer :: fail,it
      logical :: ok,raus=.false.
      character(len=2) :: bem
*-----------------------------------------------------------------------
      VIETA1(pp,qq) = qq/(-pp/2.d0-DSQRT(pp**2/4.d0-qq))
      VIETA2(pp,qq) = qq/(-pp/2.d0+DSQRT(pp**2/4.d0-qq))
      VIETA3(pp,qq) = -pp/2.d0+DSQRT(pp**2/4.d0-qq)
      VIETA4(pp,qq) = -pp/2.d0-DSQRT(pp**2/4.d0-qq)
*-----------------------------------------------------------------------
     
      !y0 = LOG(LL(0))
      !y1 = LOG(LL(1))
      !y2 = LOG(LL(2))
      !y3 = LOG(LL(3))
      !cc = 0.25*(y0-y1-y2+y3)
      !bb = (y3-y0)/3.0-3.0*cc
      !LL(4) = EXP(y0 + 4.0*bb + 16.0*cc)

      !alpha = sqrt(LL(3)*LL(0)/(LL(2)*LL(1)))
      !LL(4) = LL(3)*(LL(1)/LL(0)*alpha**3)   

      aa = L1**2-L0*L2
      bb = L3*L0-L1*L2
      cc = L2**2-L1*L3
      ok = (aa.ne.0.d0)
      bem = ''

      if (ok) then
        pp = bb/aa
        qq = cc/aa
        ok = (pp**2/4.d0.ge.qq)
        fail = 1
      endif

      if (ok) then
        if (pp>0.d0) then
          a1 = VIETA1(pp,qq)
        else  
          a1 = VIETA3(pp,qq)
        endif  
        if (a1<0.d0) then
          if (pp>0.d0) then
            a1 = VIETA4(pp,qq)
          else
            a1 = VIETA2(pp,qq)
          endif  
        endif     
        ok = (a1>0.d0.and.(L1.ne.a1*L0))
        fail = 2
      endif

      if (ok) then
        a2 = (L2-a1*L1)/(L1-a1*L0)
        if ((a2<0.2*a1).or.(a1==a2)) then
          !print*,a1,a2 
          do it=1,10
            a2 = 0.2*a1
            !aa = L1-a2*L0
            !bb = a2*(L1-a2*L0)
            !cc = a2**2*L1-L3
            aa = L2-a2**2*L0
            bb = a2*L2-L3
            cc = a2*(a2*L2-L3)
            pp = bb/aa
            qq = cc/aa
            ok = (pp**2/4.d0.ge.qq)
            if (ok) then
              a1 = VIETA1(pp,qq)
            endif  
            if ((.not.ok).or.(a1<5.0*a2)) exit
          enddo  
          ok = ok.and.(a1>0.2*a2)
          fail = 3
          bem = ' 3'
        endif  
      endif

      if (ok) then
        n2 = (L2-a1**2*L0)/(a2**2-a1**2)
        n1 = L0-n2
        ok = (n2>0.d0.and.n1>0.d0)
        fail = 4
      endif

      if (.not.ok) then
        n1 = L0
        n2 = 0.d0
        !a1 = MAX(L1/L0,(L3/L0)**(1.d0/3.d0))
        a1 = (L3/L0)**(1.d0/3.d0)
        a2 = a1
        bem = ' 4'
      endif  

      !if (.not.ok.or.raus) then
      !  print*,ip,ok,fail
      !  print*,n1,n2,a1,a2
      !  print*,L0,n1+n2
      !  print*,L1,n1*a1+n2*a2
      !  print*,L2,n1*a1**2+n2*a2**2
      !  print*,L3,n1*a1**3+n2*a2**3
      !endif  

      L1 = n1*a1+n2*a2
      L2 = n1*a1**2+n2*a2**2
      L3 = n1*a1**3+n2*a2**3
      CLOSURE = n1*a1**4 + n2*a2**4
      if (verbose>0) then
        print'(I4,2(1pE10.2),2(1pE12.4),A2)',ip,n1,n2,a1/mic,a2/mic,bem
      endif  

      end
      
