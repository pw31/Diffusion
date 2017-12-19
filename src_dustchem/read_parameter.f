************************************************************************
      subroutine READ_PARAMETER
************************************************************************
      use NATURE,ONLY: bar
      use PARAMETERS,ONLY: elements_select,model_name,dustchem_file,
     >                     struc_file,
     >                     logg,Teff,vzconst,pconst,beta,Hp,
     >                     pmin,pmax,Nl,Vl,Tfast,tsim,verbose,
     >                     bc_low,bc_high,init,implicit,tindep,
     >                     influx,outflux,inrate,outrate,vin,vout,
     >                     abund_pick,Nout,outtime,tfac,dust_diffuse
      use CHEMISTRY,ONLY: NewBackIt,NewFullIt,NewBackFac,dispol_file
      use GRID,ONLY: Npoints
      implicit none
      integer :: i,iarg,iline
      character(len=200) :: ParamFile,line

      !-------------------------
      ! ***  default values  ***
      !-------------------------
      model_name = 'output'
      struc_file = '2Drift_1800.data'
      dustchem_file = 'DustChem.dat'
      dispol_file(1) = 'dispol_BarklemCollet.dat'
      dispol_file(2) = 'dispol_StockKitzmann_withoutTsuji.dat'
      dispol_file(3) = 'dispol_WoitkeRefit.dat'
      dispol_file(4) = ''
      elements_select= 'H He C N O Na Mg Si Fe Al Ca Ti S Cl K el'
      abund_pick = 3
      Tfast      = 1000.d0
      NewFullIt  = .true.
      NewBackIt  = 5
      NewBackFac = 1.E+2
      tsim = 300.0          ! 5 minutes
      pmin       = 1.d-8*bar
      pmax       = 100.d0*bar
      beta       = 1.5
      Nl         = 1000.0
      Vl         = 3.d-23*Nl
      Npoints    = 100
      pconst     = 0.1*bar  ! 100 mbar
      vzconst    = 10.0     ! 10 cm/s
      Hp         = 1.0
      bc_low     = 1        ! fixed conc.
      bc_high    = 1        ! fixed conc.
      init       = 1        ! start with x=1 everywhere
      influx     = 0.0      
      outflux    = 0.0
      inrate     = 0.0
      outrate    = 0.0
      vin        = 1.0
      vout       = 1.0
      implicit   = .false.
      tindep     = .false.
      dust_diffuse = .false. 
      Nout       = 8
      outtime(1:8) = (/0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5/)
      tfac = 10.0
      verbose = 1

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
      iline = 0
      do 
        read(1,'(A200)',end=100) line
        if (line(1:1).eq.'#') cycle           ! ignore comment lines
        if (len(trim(line))==0) cycle         ! ignore comment lines
        iline = iline+1
        print*,trim(line)
        if (iline.eq.1) then
          elements_select = ' '//trim(line)//' '  
        else if (index(line,"! abund_pick")>0) then   
          read(line,*) abund_pick
        else if (index(line,"! Tfast")>0) then   
          read(line,*) Tfast
        else if (index(line,"! NewBackFac")>0) then 
          read(line,*) NewBackFac
        else if (index(line,"! NewBackIt")>0) then 
          read(line,*) NewBackIt
        else if (index(line,"! NewFullIt")>0) then 
          read(line,*) NewFullIt
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
        else if (index(line,"! init")>0) then   
          read(line,*) init
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
        else if (index(line,"! tindep")>0) then   
          read(line,*) tindep
        else if (index(line,"! Nout")>0) then   
          read(line,*) Nout
          read(1,*) outtime(1:Nout)           
        else if (index(line,"! tfac")>0) then   
          read(line,*) tfac
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
        else if (index(line,"! verbose")>0) then   
          read(line,*) verbose
        else if (index(line,"! dust_diffuse")>0) then   
          read(line,*) dust_diffuse
        else if (index(line,"! Nl")>0) then   
          read(line,*) Nl
          Vl = 3.d-23*Nl
        else
          print*,"*** syntax error in "//trim(ParamFile)//":"
          print*,trim(line)
          stop
        endif  
      enddo  
 100  continue
      print*,"... output directory will be "//trim(model_name)
      end
