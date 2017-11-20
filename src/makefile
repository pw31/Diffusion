F90          = ifort
F90LINKER    = ifort
#F90          = gfortran
#F90LINKER    = gfortran

# for IFORT-Compiler
#===================
DEFS      = -D IFORT -diag-disable 8290,8291 
#FFLAGS   = -r8 -g -traceback -fpe0 -check all,noarg_temp_created -warn all
FFLAGS    = -r8 -g -traceback -O3 -xHOST -prec-div -fp-model source
INCLUDES  = 
LFLAGS    = $(FFLAGS)

# for gfortran-Compiler
#======================
#DEFS      = 
#FFLAGS    = -fdefault-real-8 -fdefault-double-8 -fbounds-check -fbacktrace -g 
#FFLAGS    = -fdefault-real-8 -fdefault-double-8 -g -O5
#INCLUDES  = 
#LFLAGS    = $(FFLAGS)

EXECS    = diffuse

OBJECTS  = datamod.o \
           init_nature.o \
	   read_parameter.o \
	   read_structure.o \
           init_grid.o \
           dcoeffs.o \
           initial_conditions.o \
           diffusion.o \
	   gauss8.o \
	   linpack_q.o \
           main.o

.SUFFIXES : .o .f .F .f90

default: $(OBJECTS)
	 $(F90LINKER) -o $(EXECS) $(OBJECTS) $(LFLAGS) 

clean:
	rm -f *.o 

realclean:
	rm -f *.o *~ *.mod *__genmod.f90 $(EXECS)

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
.F.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

