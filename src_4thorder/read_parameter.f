************************************************************************
      subroutine READ_PARAMETER
************************************************************************
      use NATURE,ONLY: bar
      use PARAMETERS,ONLY: logg,Teff,vzconst,pconst,beta,Hp,pmin,pmax,
     >                     bc_low,bc_high,init,implicit,tindep,simple,
     >                     influx,outflux,inrate,outrate,vin,vout,
     >                     Nout,outtime,tfac
      use READMODEL,ONLY: struc_file
      use GRID,ONLY: Npoints
      implicit none
      integer :: i,iarg,iline
      character(len=200) :: ParamFile,line

      !-------------------------
      ! ***  default values  ***
      !-------------------------
      struc_file = '2Drift_1800.data'
      pmin = 1.d-8*bar
      pmax = 100.0*bar
      Npoints = 100
      pconst  = 0.1*bar  ! 100 mbar
      vzconst = 10.0     ! 10 cm/s
      beta = 1.5
      Hp = 1.0
      bc_low  = 1        ! fixed conc.
      bc_high = 1        ! fixed conc.
      init    = 1        ! start with x=1 everywhere
      influx  = 0.0      
      outflux = 0.0
      inrate  = 0.0
      outrate = 0.0
      vin     = 1.0
      vout    = 1.0
      simple   = .true.
      implicit = .false.
      tindep   = .false.
      Nout = 8
      outtime(1:8) = (/0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5/)
      tfac = 10.0

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
        if (index(line,"! Npoints")>0) then   
          read(line,*) Npoints
        else if (index(line,"! struc_file")>0) then
          i = index(line,' ')
          struc_file = line(1:i-1)
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
        else if (index(line,"! simple")>0) then   
          read(line,*) simple
        else if (index(line,"! Nout")>0) then   
          read(line,*) Nout
          read(1,*) outtime(1:Nout)           
        else if (index(line,"! tfac")>0) then   
          read(line,*) tfac
        else
          print*,"*** syntax error in "//trim(ParamFile)//":"
          print*,trim(line)
          stop
        endif  
      enddo  
 100  continue
      end
