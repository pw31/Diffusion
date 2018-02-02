!----------------------------------------------------------------------
      SUBROUTINE INIT_NATURE
!----------------------------------------------------------------------
      use NATURE,ONLY: cl,hplanck,bk,elad,grav,NA,pi,sig_SB,Rgas,
     >                 Msun,Mearth,MJup,Lsun,Rsun,amu,mel,yr,km,AU,
     >                 Ang,nm,mic,pc,eV,Ws,Wm2Hz,Jansky,Mbarn,bar,atm
      implicit none
      real :: cPl1, cPl2

!     -----------------------------------------------------
!     ***        fundamental nature constants           *** 
!     ***              [cgs] (from FLASH)               *** 
!     -----------------------------------------------------
      cl      = 2.9979245800000E+10 ! speed of light
      hplanck = 6.6260755400000E-27 ! Planck's constant
      bk      = 1.3806581200000E-16 ! Boltzmann's constant
      elad    = 4.8032068150000E-10 ! electron charge
      grav    = 6.6725985000000E-08 ! gravitational constant
      NA      = 6.0221417900000E+23 ! Avogadro constant
      Rgas    = 8.3144598000000E+00 ! ideal gas constant [J/mol/K]
      pi      = ACOS(-1.d0)

!     ---------------------------
!     ***  derived constants  ***
!     ---------------------------
      cPl1   = 2.0*hplanck*cl**2
      cPl2   = hplanck*cl/bk    
      sig_SB = cPl1/cPl2**4*pi**5/15.0

!     ---------------------------
!     ***  other units [cgs]  ***
!     ---------------------------
      Msun   = 1.988922500E+33     ! solar mass
      Mearth = 5.974200000E+27     ! mass of Earth
      Mjup   = 1.898600000E+30     ! mass of Jupiter 
      Lsun   = 3.846000000E+33     ! solar luminosity
                                   ! (NASA solar physics division)
      Rsun   = 6.959900000E+10     ! solar radius     
                                   ! (Brown & Christensen-Dalsgaard 1998)
      amu    = 1.660531000E-24     ! atomar mass unit
      mel    = 9.109389754E-28     ! electron mass
      yr     = 3.155760000E+7      ! 1 year
      km     = 1.E+5               ! 1 kilometer
      AU     = 1.495978700E+13     ! 1 astronomical unit
      Ang    = 1.E-8               ! 1 Angstroem
      nm     = 1.E-7               ! 1 nanometer
      mic    = 1.E-4               ! 1 micron
      pc     = 3.085680250E+18     ! 1 parsec
      eV     = 1.602176462E-12     ! 1 electron Volt
      Ws     = 1.E+7               ! 1 Watt*sec = 1 Joule 
      Wm2Hz  = 1.E+3               ! 1 W/m^2/Hz in [erg/cm^2/s/Hz]
      Jansky = 1.E-23              ! 1 Jy in [erg/cm^2/s/Hz]
      Mbarn  = 1.E-18              ! 1 MegaBarn in [cm^2] 
      bar    = 1.E+6               ! 1 bar in [dyn/cm2]
      atm    = 1.013D+6            ! standard atmosphere pressure

      RETURN
      end


