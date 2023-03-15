# v. 04.10.2009 for stallo - copied from phys version of cluster v.3.3

# Ver. 1.0. 27/8-10 - Based on Makefile for Einar Stiansen.

# If using Intel compilers, uncomment a line at the top of the igrf.f90 file
FC 	= gfortran  # mpif90, ifort
LD      = $(FC)
LDFLAGS = #-ipo
FFLAGS 	= -Wall -O3 -fPIC -fcheck=all -fbounds-check -pedantic -Wextra #-ipo -check bounds -check pointer -check uninit
SFLAGS  = -static #-libgfortran -static-libgcc (AGB, static should work?)
STATIC  = 0

ifneq ($(STATIC),0)
	LDFLAGS += $(SFLAGS)
endif

#MKL_INCLUDE = -I/global/apps/intel/mkl/10.1.0.015/include/
LIBS 	= #-L/global/apps/intel/mkl/10.1.0.015/ -shared-intel -Wl,--start-group /global/apps/intel/mkl/10.1.0.015/lib/em64t/libmkl_intel_lp6\4.a /global/apps/intel/mkl/10.1.0.015/lib/em64t/libmkl_core.a

PROG	   = apextest

OBJS	   = igrf.o magfld.o checkapexsh.o apexsh.o apex.o makeapexsh.o

$(PROG)	    : $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(PROG) $(LIBS)

apexsh.o      : apexsh.f90
	$(FC) -c $(FFLAGS) $<

igrf.o      : igrf.f90
	$(FC) -c $(FFLAGS) $<

checkapexsh.o        : checkapexsh.f90
	$(FC) -c $(FFLAGS) $<


apex.o      : apex.f90
	$(FC) -c $(FFLAGS) $<

magfld.o      : magfld.f90
	$(FC) -c $(FFLAGS) $<

makeapexsh.o      : makeapexsh.f90 apexsh.o
	$(FC) -c $(FFLAGS) $<



clean:
	rm -f *.o *.mod *.dat apextest
