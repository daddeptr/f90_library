# Config file for linking against packages in
# HPCPorts.  Should work on all NERSC machines
# after loading the cmb module from HPCPorts.

# =============== Installation ===================

# Final install location.  This defaults to a
# subdirectory of the top-level directory called
# install_$COMMANDER.
#export INSTALL := /usr/local


# =============== System Tools ===================

# The defaults should be fine for any normal
# GNU/Linux operating system...

# Here you can override the make version to use
#export MAKE := gmake

# How to make a tar.gz file
#export CTAR := tar czvf

# The archive command
#export AR := ar

# Archive options to create a library
#export ARFLAGS := crv

# The ranlib command
#export RANLIB := ranlib


# ================= Compilers ====================

# Serial Fortran 90 compiler.  Default is "f90".
F90=$(HPCP_FC)


# ============== Language Mixing =================

 MPFCLIBS=$(HPCP_MPIFCLIBS) $(HPCP_FCLIBS)


# =========== Compiler Optimizations =============

# FIXME! the line length option is gfortran specific, but I 
# gave up hacking the code to make it conform to the standard
# line length limit.

# Fortran 90 compiler flags
F90FLAGS=$(HPCP_FCFLAGS) $(HPCP_OMPFLAGS) -ffree-line-length-none -fbacktrace -g

# ============== Fortran Features ================

# Set this variable to 1 if the fortran compiler
# produces module files with capitalization (*.MOD)
# instead of lowercase (*.mod).
#export FORTRAN_UPPER := 1


# =============== CFITSIO Linking ================

# The include and linking commands for cfitsio.
CFITSIO_LINK=$(cfitsio_LIBS_FC)


# =============== LAPACK Linking =================

LAPACK_LINK=$(lapack_LIBS_FC) $(blas_LIBS_FC)


# ================ Healpix linking ================

HEALPIX_INCLUDE=-I$(healpix_PREFIX)/include
HEALPIX_LINK=$(healpix_LIBS_FC)

all: main

main: 
	$(F90) $(F90FLAGS) $(CFITSIO_LINK) $(HEALPIX_INCLUDE) $(HEALPIX_LINK) alm_mcmc.f90 -o alm_mcmc.x
