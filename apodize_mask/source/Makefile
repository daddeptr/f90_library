Healpix_LIB = /global/scratch2/sd/dpietrob/Software/Healpix_3.00/lib_ifort
### Healpix_LIB = /Users/dpietrob/work/cmbcodes/Healpix_3.10/lib_gfortran

Healpix_INCLUDE = /global/scratch2/sd/dpietrob/Software/Healpix_3.00/include_ifort
### Healpix_INCLUDE = /Users/dpietrob/work/cmbcodes/Healpix_3.10/include_gfortran

CFITSIO = /global/scratch2/sd/dpietrob/Software/cfitsio_ifort
### CFITSIO = /opt/local/lib

F90C    = ifort
FFLAGS  = -mkl -openmp -O4 -xHost -no-prec-div -fpp -DMPI -CB $(CHECKFLAGS)
CHECKFLAGS = -check all -warn -CB -traceback -g

### gfortran
### F90C    = gfortran
### CHECKFLAGS = -fbounds-check -g -fbacktrace -pedantic -Wall -pedantic-errors -fcheck-array-temporaries -finit-local-zero
### FFLAGS  = -O3 -I$(F90_INCDIR) -DGFORTRAN -fno-second-underscore -fopenmp -ffree-line-length-none $(CHECKFLAGS)
### FFLAGS  = -O3 -DGFORTRAN -fno-second-underscore -m64 -fopenmp -ffree-line-length-none #$(CHECKFLAGS)

LIBRARIES = -L$(Healpix_LIB) -L$(CFITSIO) -lhealpix -lcfitsio

LAPACKL =

NAME = apodize_mask.x

EXECDIR = ./

F90FLAGS = -I$(Healpix_INCLUDE) $(FFLAGS) 

OBJFILES = apodize_mask.o

default: main

.f.o:
	f77 $(F90FLAGS) -c $<

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $*.F90


main: $(OBJFILES) 	
	mkdir -p $(EXECDIR)
	$(F90C) $(F90FLAGS) $(OBJFILES) $(LIBRARIES) -o $(EXECDIR)/$(NAME)

all: main

clean:
	rm -f *.o *.mod *.d *.pc *.obj core
tidy: clean
	rm $(EXECDIR)/$(NAME)
