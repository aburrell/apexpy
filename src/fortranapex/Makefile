# v. 04.10.2009 for stallo - copied from phys version of cluster v.3.3

# Ver. 1.0. 27/8-10 - Based on Makefile for Einar Stiansen.

FC 	= gfortran
#FC 	= mpif90
LD      = $(FC)
LDFLAGS = #-ipo
FFLAGS 	= -O3 -fPIC #-ipo -check bounds -check pointer -check uninit

#MKL_INCLUDE = -I/global/apps/intel/mkl/10.1.0.015/include/
LIBS 	= #-L/global/apps/intel/mkl/10.1.0.015/ -shared-intel -Wl,--start-group /global/apps/intel/mkl/10.1.0.015/lib/em64t/libmkl_intel_lp6\4.a /global/apps/intel/mkl/10.1.0.015/lib/em64t/libmkl_core.a

PROG	   = apextest

OBJS	   = checkapexsh.o apexsh.o makeapexsh.o apex.o magfld.o

$(PROG)	    : $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(PROG) $(LIBS)

apexsh.o      : apexsh.f90  
	$(FC) -c $(FFLAGS) $<

checkapexsh.o        : checkapexsh.f90 
	$(FC) -c $(FFLAGS) $<


apex.o      : apex.f  
	$(FC) -c $(FFLAGS) $<

magfld.o      : magfld.f
	$(FC) -c $(FFLAGS) $<

makeapexsh.o      : makeapexsh.f90 apexsh.o  
	$(FC) -c $(FFLAGS) $<



clean:
	rm -f *.o *.mod *.dat apextest
