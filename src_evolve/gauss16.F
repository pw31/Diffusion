**********************************************************************
           SUBROUTINE GAUSS16 (Ndim,N,a,x,b)
**********************************************************************
*****                                                            *****
*****   Diese Routine loesst ein lineares Gleichungssystem       *****
*****   der Form    (( a )) * ( x ) = ( b )     nach x auf.      *****
*****   Der Algorithmus funktioniert, indem die Matrix a         *****
*****   auf Dreiecksform gebracht wird.                          *****
*****                                                            *****
*****   EINGABE:  n = Dimension der Vektoren, der Matrix         *****
*****             a = (N x N)-Matrix                             *****
*****             b = (N)-Vektor                                 *****
*****   AUSGABE:  x = (N)-Vektor                                 *****
*****                                                            *****
**********************************************************************
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer :: Ndim,N,i,j,k,kmax
      real(kind=qp) :: a(Ndim,Ndim),x(Ndim),b(Ndim),c,amax

      do i=1,N-1
        !------------------------------------------
        !***  MAX-Zeilentausch der i-ten Zeile  ***      
        !------------------------------------------
        kmax = i
        amax = ABS(a(i,i))
        do k=i+1,N
          if (ABS(a(k,i))>amax) then
            kmax = k
            amax = ABS(a(k,i))
          endif
        enddo
        if (kmax.ne.i) then
          do j=1,N
            c         = a(i,j)
            a(i,j)    = a(kmax,j)
            a(kmax,j) = c 
          enddo  
          c       = b(i)
          b(i)    = b(kmax)
          b(kmax) = c 
        endif
        !---------------------------------
        !***  bringe auf Dreiecksform  ***
        !---------------------------------
        do k=i+1,N
          c = a(k,i) / a(i,i)
          a(k,i) = 0.Q0
          do j=i+1,N
            a(k,j) = a(k,j) - c*a(i,j)
          enddo  
          b(k) = b(k) - c * b(i)
        enddo  
      enddo  

      !--------------------------
      !***  loese nach x auf  ***
      !--------------------------
      do i=N,1,-1
        c = 0.Q0
        if (i.lt.N) then
          do j=i+1,N
            c = c + a(i,j) * x(j)
          enddo
        end if
        x(i) = (b(i) - c) / a(i,i)
      enddo
  
      RETURN
      end
