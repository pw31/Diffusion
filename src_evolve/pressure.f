      real*8 FUNCTION PRESSURE(Tg)
      use NATURE,ONLY: bk
      use ELEMENTS,ONLY: NELEM
      use CHEMISTRY,ONLY: NMOLE
      use EXCHANGE,ONLY: nmol,nat,nel
      implicit none
      real*8,intent(in) :: Tg
      integer :: j
      real*8 :: kT,nges

      kT = bk*Tg
      nges = nel
      do j=1,NELEM
        nges = nges + nat(j)
      enddo
      do j=1,NMOLE
        nges = nges + nmol(j)    ! includes the ions
      enddo
      PRESSURE = nges*kT

      end
