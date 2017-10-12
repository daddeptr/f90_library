
#CAMB Makefile
#Edit for your compiler
#Note there are many ifc versions, some of which behave oddly

#Intel 8 , -openmp toggles mutli-processor:
F90C     = ifort
#FFLAGS = -O2 -Vaxlib -W0 -WB -fpp -ip -i-static -openmp
FFLAGS = -check all -check bounds -traceback #-i-statis
FFlAGS = -O3 -fast #-O2 -Vaxlib -openmp -fpp2 -lmpi -lscs  -ldl -DMPI

EQUATIONS     = equations3
POWERSPECTRUM = power_tilt
#Module doing non-linear scaling
NONLINEAR     = halofit

#Driver program
DRIVER        = inidriver2.F90
#DRIVER        = sigma8.f90
#DRIVER        = tester.f90

#Settings for building camb_fits
#Location of FITSIO and name of library
#FITSDIR       = /home/cpac/cpac-tools/lib
FITSLIB       = cfitsio
#Location of HEALPIX for building camb_fits
HEALPIXDIR    = /usr/local/phys_progr/Healpix_2.01/lib

CAMBLIB       = ../libcamb.a
CAMBNAME      = ../qexp2Camb

#Shouldn't need to change anything else...

F90FLAGS      = $(FFLAGS)
HEALPIXLD     = -L$(HEALPIXDIR)/lib -lhealpix -L$(FITSDIR) -l$(FITSLIB)
FC            = $(F90C)

CAMBOBJ       = subroutines.o inifile.o $(POWERSPECTRUM).o recfast.o modules2.o \
	bessels.o mysub.o $(EQUATIONS).o $(NONLINEAR).o lensing.o cmbmain.o camb2.o 

default: camb

all: camb $(CAMBLIB)

camb: $(CAMBOBJ) $(DRIVER)
	$(F90C) $(F90FLAGS) $(CAMBOBJ) $(DRIVER) -o $(CAMBNAME)

$(CAMBLIB): $(CAMBOBJ)
	ar -r $@ $?

camb_fits: writefits.f90 $(CAMBOBJ) $(DRIVER)
	$(F90C) $(F90FLAGS) -I$(HEALPIXDIR)/include $(CAMBOBJ) writefits.f90 $(DRIVER) $(HEALPIXLD) -DWRITE_FITS -o $@

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90

clean:
	-rm -f *.o *.a *.d core *.mod $(CAMBNAME)

tidy:
	-rm -f *.o *.a *.d core *.mod $(CAMBNAME) $(CAMBLIB)
