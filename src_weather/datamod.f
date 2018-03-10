************************************************************************
      module NATURE
************************************************************************
      real :: cl,hplanck,bk,elad,grav,NA,pi,sig_SB,Rgas,
     >        Msun,Mearth,MJup,Lsun,Rsun,amu,mel,yr,km,AU,
     >        Ang,nm,mic,pc,eV,Ws,Wm2Hz,Jansky,Mbarn,bar,atm
      end

************************************************************************
      module PARAMETERS
************************************************************************
      character(len=200) :: elements_select,model_name,dustchem_file
      character(len=200) :: struc_file
      real :: logg,Teff,vzconst,pconst,beta,Hp,pmin,pmax,Nl,Vl
      real :: influx,outflux,inrate,outrate,vin,vout
      real :: Tfast,tfac,outtime(50),tsim,dtfix,precision
      integer :: bc_low,bc_high,init,Nout,abund_pick,evap_model,verbose
      logical :: implicit,tindep,dust_diffuse
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
      real,allocatable :: BB(:,:)               ! implicit diffusion matrix
      real,allocatable :: xlower(:),xupper(:)   ! boudary values
      real :: dt_diff_ex,dt_diff_im
      end

************************************************************************
      module STRUCT
************************************************************************
      real,allocatable :: Diff(:)     ! diffusion coefficient [cm2/s]
      real,allocatable :: rho(:)      ! mass density [g/cm3]
      real,allocatable :: nHtot(:)    ! total H nuclei density [1/cm3]
      real,allocatable :: Temp(:)     ! temperature [K]
      real,allocatable :: press(:)    ! pressure [dyn/cm2]
      real,allocatable :: mu(:)       ! mean molecular weight [g]
      real,allocatable :: nHeps(:,:)  ! element abundance [cm-3]
      real,allocatable :: rhoLj(:,:)  ! dust moments [cm^(j-3)]
      real,allocatable :: rhoL3(:,:)  ! dust volumes [cm^3/cm^3]
      real,allocatable :: mols(:,:)   ! molecular particle densities
      real,allocatable :: atms(:,:)   ! atomic particle densities
      real,allocatable :: elec(:)     ! electron particle densities
      end

************************************************************************
      module ELEMENTS
************************************************************************
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,parameter :: NELEM=41           ! number of elements (up to Zr + W)
      integer :: NEPS                         ! number of affected elements
      character(len=2) :: elnam(NELEM)        ! names of elements
      real(kind=qp) :: eps0(NELEM)            ! element abundances
      integer :: elnr(NELEM),elcode(NELEM)    ! element cross-indices
      real*8 :: mass(NELEM)                   ! element masses
      real*8 :: muH                           ! rho/n<H>
      end

************************************************************************
      module DUST_DATA
************************************************************************
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,parameter :: NDUSTmax=200       ! max number of condensed species
      integer :: NDUST                        ! number of condensed species      
      character(len=20) :: dust_nam(NDUSTmax) ! names of dust species
      real*8  :: dust_rho(NDUSTmax)           ! dust material densities
      real*8  :: dust_mass(NDUSTmax)          ! dust monomer volume
      real*8  :: dust_vol(NDUSTmax)           ! dust monomer volume
      real*8  :: Tmelt(NDUSTmax)              ! melting points
      real*8  :: Tcorr(NDUSTmax)
      logical :: is_liquid(NDUSTmax)
      integer :: dust_nel(NDUSTmax)           ! no of elements in dust
      integer :: dust_el(NDUSTmax,8)          ! indices of elements
      integer :: dust_nu(NDUSTmax,8)          ! stoichiometric coeffs      
      integer :: fit(NDUSTmax)                ! fit-formular identifier
      real*8  :: cfit(NDUSTmax,0:4)           ! pvap fit coefficients
      end

************************************************************************
      module CHEMISTRY
************************************************************************
      use ELEMENTS,ONLY: NELEM
      character(len=200) :: dispol_file(4)
      logical :: NewFullIt
      real*8  :: NewBackFac
      integer :: NewBackIt,NewFastLevel
      integer :: NMOLdim         ! max number of molecules
      integer :: NMOLE           ! number of molecules found
      integer :: NELM            ! number of elements found
      integer :: el=0,H=0,He=0,Li=0,Be=0,B=0,C=0,N=0,O=0,F=0,Ne=0
      integer :: Na=0,Mg=0,Al=0,Si=0,P=0,S=0,Cl=0,Ar=0,K=0,Ca=0
      integer :: Sc=0,Ti=0,V=0,Cr=0,Mn=0,Fe=0,Co=0,Ni=0,Cu=0,Zn=0
      integer :: Ga=0,Ge=0,As=0,Se=0,Br=0,Kr=0,Rb=0,Sr=0
      integer :: Y=0,Zr=0,W=0
      logical :: charge
      character(len=2) :: catm(NELEM)           ! names of elements
      character(len=20),allocatable :: cmol(:)  ! names of molecules
      real,allocatable :: molmass(:)            ! masses of molecules
      integer :: elnum(NELEM)                   ! indices of found elements
      integer :: elion(NELEM)                   ! indices of ions
      integer,allocatable :: fit(:)             ! fit-formular identifier
      integer,allocatable :: natom(:)           ! no of atoms in molecule    
      integer,allocatable :: source(:)          ! no of source file
      integer,allocatable :: m_kind(:,:)        ! index of elements
      integer,allocatable :: m_anz(:,:)         ! stoichiometric coeffs
      real*8,allocatable  :: a(:,:)             ! kp fit-coeffs
      real*8,allocatable  :: error(:)           ! kp fit errors
      real*8 :: th1,th2,th3,th4,TT1,TT2,TT3     
!$omp threadprivate(th1,th2,th3,th4,TT1,TT2,TT3) 
      end

************************************************************************
      module NUCLEATION
************************************************************************
      integer :: NNUC
      integer,allocatable :: nuc(:)
      character(len=10),allocatable :: nuc_nam(:)
      end

************************************************************************
      module SPECIES
************************************************************************
      integer :: NSPECIES
      character(len=20),allocatable :: spnam(:)
      logical,allocatable :: keysp(:)
      integer,allocatable :: spnr(:),spzahl(:)
      real,allocatable    :: spmass(:)
      integer,allocatable :: spelem(:,:),spstoich(:,:)
      end
      
************************************************************************
      module REACTIONS
************************************************************************
      integer :: NREAC
      integer,allocatable :: neduct(:),nprod(:)
      integer,allocatable :: reac_sp(:,:,:),reac_nu(:,:,:)
      real,allocatable :: stick(:)
      end

************************************************************************
      module EXCHANGE
************************************************************************
      use CHEMISTRY,ONLY: NMOLE
      use ELEMENTS,ONLY: NELEM
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: nel,nat(NELEM),nion(NELEM)
      real(kind=qp),allocatable :: nmol(:)
      real*8,allocatable :: Jst(:),chi(:)
      integer :: HII,HeII,CII,NII,OII,NaII,MgII,LiII,ClII
      integer :: AlII,KII,TiII,SII,SiII,FeII,CaII,ipoint
      integer,parameter :: H=1,He=2,Li=3,Be=4,B=5,C=6,N=7,O=8,F=9
      integer,parameter :: Ne=10,Na=11,Mg=12,Al=13,Si=14,P=15,S=16
      integer,parameter :: Cl=17,Ar=18,K=19,Ca=20,Sc=21,Ti=22
      integer,parameter :: V=23,Cr=24,Mn=25,Fe=26,Co=27,Ni=28
      integer,parameter :: Cu=29,Zn=30,Ga=31,Ge=32,As=33,Se=34
      integer,parameter :: Br=35,Kr=36,Rb=37,Sr=38,Y=39,Zr=40,W=41
      integer*8 :: chemcall=0,chemiter=0,itransform=0,ieqcond=0
      integer*8 :: Fcall=0,Jcall=0
      logical,allocatable :: inactive(:)
!$omp threadprivate(nel,nat,nion,nmol,Jst,chi,ipoint,inactive)
      end
