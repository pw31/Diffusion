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
      real :: logg,Teff,vzconst,pconst,beta,Hp,pmin,pmax,tnull,tend
      real :: influx,outflux,inrate,outrate,vin,vout
      real :: tfac,outtime(50)
      integer :: bc_low,bc_high,init,Nout
      logical :: simple,implicit,tindep
      end

************************************************************************
      module READMODEL
************************************************************************
      character(len=200) :: struc_file
      integer :: Nlayers
      real,dimension(1000) :: Rlay,Tlay,play,rholay,glay
      real,dimension(1000) :: zlay,mulay,vconvlay,Difflay
      end

************************************************************************
      module GRID
************************************************************************
      integer :: Npoints                        ! number of gridpoints
      real,allocatable :: zz(:)                 ! vertical gridpoints
      real,allocatable :: d1l2(:),d1l1(:),d1m(:),d1r1(:),d1r2(:)  ! first derivatives
      real,allocatable :: d2l2(:),d2l1(:),d2m(:),d2r1(:),d2r2(:)  ! second derivatives
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
