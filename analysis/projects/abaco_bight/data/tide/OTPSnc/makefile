ARCH = $(shell uname -s)
ifeq ($(ARCH),Linux)
 FC = gfortran
 NCLIB = /usr/local/lib
 NCINCLUDE = /usr/local/include
 NCLIBS= -lnetcdf -lnetcdff
endif

predict_tide: predict_tide.f90 subs.f90 constit.h
	$(FC) -o predict_tide predict_tide.f90 subs.f90 -L$(NCLIB) $(NCLIBS) -I$(NCINCLUDE) 
	#rm *.o
extract_HC:  extract_HC.f90 subs.f90
	$(FC) -o  extract_HC extract_HC.f90 subs.f90 -L$(NCLIB) $(NCLIBS) -I$(NCINCLUDE)
	#rm *.o
