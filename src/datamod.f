************************************************************************
      module NATURE
************************************************************************
      real :: cl,hplanck,bk,elad,grav,NA,pi,sig_SB,Rgas,
     >        Msun,Mearth,MJup,Lsun,Rsun,amu,mel,yr,km,AU,
     >        Ang,nm,mic,pc,eV,Ws,Wm2Hz,Jansky,Mbarn,bar
      end

************************************************************************
      module PARAMETERS
************************************************************************
      real :: logg,Teff,vzconst,pconst,beta,Hp,hmin,hmax,tnull,tend
      real :: influx,outflux,inrate,outrate,vin,vout
      real :: tfac,outtime(50)
      integer :: bc_low,bc_high,init,Nout
      logical :: simple,implicit,tindep
      end

************************************************************************
      module READMODEL
************************************************************************
      integer :: Nlayers
      real,dimension(1000) :: Rlay,Tlay,play,rholay,glay
      real,dimension(1000) :: zlay,mulay,vconvlay,Difflay
      end

************************************************************************
      module GRID
************************************************************************
      integer :: Npoints                        ! number of gridpoints
      real,allocatable :: zz(:)                 ! vertical gridpoints
      real,allocatable :: d1l(:),d1m(:),d1r(:)  ! first derivatives
      real,allocatable :: d2l(:),d2m(:),d2r(:)  ! second derivatives
      real,allocatable :: xx(:)                 ! concentration
      end

************************************************************************
      module STRUCT
************************************************************************
      real,allocatable :: Diff(:)     ! diffusion coefficient [cm2/s]
      real,allocatable :: rho(:)      ! mass density [g/cm3]
      real,allocatable :: nHtot(:)    ! total H nuclei density [1/cm3]
      real,allocatable :: T(:)        ! temperature [K]
      real,allocatable :: press(:)    ! pressure [dyn/cm2]
      real,allocatable :: mu(:)       ! mean molecular weight [g]
      real,allocatable :: eps(:)      ! element abundance [-]
      end
