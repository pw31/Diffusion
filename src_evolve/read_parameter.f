************************************************************************
      subroutine READ_PARAMETER
************************************************************************
      use NATURE,ONLY: bar,yr,km
      use PARAMETERS,ONLY: elements_select,model_name,dustchem_file,
     >                     struc_file,Tfast,Tcrust,
     >                     tsim,dt_init,dt_increase,dt_max,heatrate,
     >                     logg,vzconst,pconst,beta,Hp,Rplanet,
     >                     gas_kind,crust_kind,crust_thickness,
     >                     pmin,pmax, bc_low,bc_high,implicit,
     >                     influx,outflux,inrate,outrate,vin,vout,
     >                     useDatabase,verbose,immediateEnd
      use CHEMISTRY,ONLY: NewBackIt,NewFullIt,NewBackFac,NewFastLevel,
     >                    dispol_file
      use GRID,ONLY: Npoints
      implicit none
      integer :: i,iarg,iline,dispol_set
      character(len=200) :: ParamFile,line

      !-------------------------
      ! ***  default values  ***
      !-------------------------
      model_name = 'output'
      struc_file = 'none'
      dustchem_file = 'DustChem.dat'
      dispol_file(1) = 'dispol_BarklemCollet.dat'
      dispol_file(2) = 'dispol_StockKitzmann_withoutTsuji.dat'
      dispol_file(3) = 'dispol_WoitkeRefit.dat'
      dispol_file(4) = ''
      elements_select= 'H He C N O Na Mg Si Fe Al Ca Ti S Cl K el'
      Tfast      = 1000.d0
      NewFullIt  = .true.
      NewBackIt  = 5
      NewBackFac = 1.E+2
      NewFastLevel = 0
      tsim       = 300.0    ! 5 minutes
      dt_init    = 1.E-3    ! 1 milli-sec
      dt_increase= 1.3      ! factor for dt increase
      dt_max     = 9.E+99   ! unlimited timestep
      heatrate   = 0.0      ! heating rate [K/yr]
      pmin       = 1.d-8*bar
      pmax       = 100.d0*bar
      beta       = 1.5
      Npoints    = 100
      pconst     = 0.1*bar    ! 100 mbar
      vzconst    = 10.0       ! 10 cm/s
      Hp         = 1.0
      Tcrust     = 1000.0     ! surface temperature [K]
      Rplanet    = 6371.0*km  ! planet radius [cm]
      logg       = 9.81*100.0 ! surface gravity [cm/s2]
      bc_low     = 1          ! fixed conc.
      bc_high    = 1          ! fixed conc.
      influx     = 0.0      
      outflux    = 0.0
      inrate     = 0.0
      outrate    = 0.0
      vin        = 1.0
      vout       = 1.0
      implicit   = .true.
      gas_kind = 1             ! 0=empty, 1=solar, 2=meteor, 3=EarthCrust
      crust_kind = 1           ! 1=solar, 2=meteor, 3=EarthCrust
      crust_thickness = 100.0  ! 1 metre
      useDatabase = .true.
      verbose = 1
      immediateEnd = .false.

      !-------------------------------------------
      ! ***  change parameters via input file  ***
      !-------------------------------------------
      iarg = iargc()
      if (iarg==0) then
        print*,"using default parameters"
        return
      endif  
      call getarg(1,ParamFile)
      print*,"reading "//trim(ParamFile)//" ..."
      open(unit=1,file=ParamFile,status='old')
      dispol_set = 0
      iline = 0
      do 
        read(1,'(A200)',end=100) line
        if (line(1:1).eq.'#') cycle           ! ignore comment lines
        if (len(trim(line))==0) cycle         ! ignore comment lines
        iline = iline+1
        print*,trim(line)
        if (iline.eq.1) then
          elements_select = ' '//trim(line)//' '  
        else if (index(line,"! Tfast")>0) then   
          read(line,*) Tfast
        else if (index(line,"! NewBackFac")>0) then 
          read(line,*) NewBackFac
        else if (index(line,"! NewBackIt")>0) then 
          read(line,*) NewBackIt
        else if (index(line,"! NewFullIt")>0) then 
          read(line,*) NewFullIt
        else if (index(line,"! NewFastLevel")>0) then 
          read(line,*) NewFastLevel
        else if (index(line,"! Npoints")>0) then   
          read(line,*) Npoints
        else if (index(line,"! pmin")>0) then   
          read(line,*) pmin
        else if (index(line,"! pmax")>0) then   
          read(line,*) pmax
        else if (index(line,"! vzconst")>0) then   
          read(line,*) vzconst
        else if (index(line,"! pconst")>0) then   
          read(line,*) pconst
          pconst = pconst*bar
        else if (index(line,"! beta")>0) then   
          read(line,*) beta
        else if (index(line,"! Hp")>0) then   
          read(line,*) Hp
        else if (index(line,"! bc_low")>0) then   
          read(line,*) bc_low
        else if (index(line,"! bc_high")>0) then   
          read(line,*) bc_high
        else if (index(line,"! gas_kind")>0) then   
          read(line,*) gas_kind
        else if (index(line,"! crust_kind")>0) then   
          read(line,*) crust_kind
        else if (index(line,"! crust_thickness")>0) then   
          read(line,*) crust_thickness
        else if (index(line,"! influx")>0) then   
          read(line,*) influx
        else if (index(line,"! outflux")>0) then   
          read(line,*) outflux
        else if (index(line,"! inrate")>0) then   
          read(line,*) inrate
        else if (index(line,"! outrate")>0) then   
          read(line,*) outrate
        else if (index(line,"! vin")>0) then   
          read(line,*) vin
        else if (index(line,"! vout")>0) then   
          read(line,*) vout
        else if (index(line,"! implicit")>0) then   
          read(line,*) implicit
        else if (index(line,"! model_name")>0) then
          i = index(line,' ')
          model_name = 'output_'//line(1:i-1)
        else if (index(line,"! dustchem_file")>0) then
          i = index(line,' ')
          dustchem_file = line(1:i-1)
        else if (index(line,"! struc_file")>0) then
          i = index(line,' ')
          struc_file = line(1:i-1)
        else if (index(line,"! tsim")>0) then   
          read(line,*) tsim
        else if (index(line,"! dt_init")>0) then   
          read(line,*) dt_init
        else if (index(line,"! dt_increase")>0) then   
          read(line,*) dt_increase
        else if (index(line,"! dt_max")>0) then   
          read(line,*) dt_max
        else if (index(line,"! heatrate")>0) then   
          read(line,*) heatrate
          heatrate = heatrate/yr
        else if (index(line,"! verbose")>0) then   
          read(line,*) verbose
        else if (index(line,"! immediateEnd")>0) then   
          read(line,*) immediateEnd
        else if (index(line,"! Tcrust")>0) then
          read(line,*) Tcrust
        else if (index(line,"! Rplanet")>0) then
          read(line,*) Rplanet
          Rplanet = Rplanet*km
        else if (index(line,"! logg")>0) then
          read(line,*) logg
        else if (index(line,"! dispol_file2")>0) then 
          i = index(line,"!")
          read(line(1:i-1),*)  dispol_file(2)
          dispol_set = 2
        else if (index(line,"! dispol_file3")>0) then 
          i = index(line,"!")
          read(line(1:i-1),*)  dispol_file(3)
          dispol_set = 3
        else if (index(line,"! dispol_file4")>0) then 
          i = index(line,"!")
          read(line(1:i-1),*)  dispol_file(4)
          dispol_set = 4
        else if (index(line,"! dispol_file")>0) then 
          i = index(line,"!")
          read(line(1:i-1),*)  dispol_file(1)
          dispol_set = 1
        else
          print*,"*** syntax error in "//trim(ParamFile)//":"
          print*,trim(line)
          stop
        endif  
      enddo  
 100  continue
      if (dispol_set>0.and.dispol_set<4) dispol_file(4)=""
      if (dispol_set>0.and.dispol_set<3) dispol_file(3)=""
      if (dispol_set>0.and.dispol_set<2) dispol_file(2)=""

      print*,"... output directory will be "//trim(model_name)
      end
