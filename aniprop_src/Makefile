FC=gfortran
CC=cc
FFLAGS=-c -g -O5 


OBJS=main.o forward_solver.o \
	 aniprop_subroutines.o \
	 rf_aniso_subroutines.o \
	 util_subroutines.o \
	 refft.o \
	 eispack.o

TARGETS=main

all : $(TARGETS) 

main : $(OBJS)
	$(FC) -o main $(OBJS) 

main.o : main.f90
	$(FC) $(FFLAGS)  main.f90 

forward_solver.o : forward_solver.f90 
	$(FC) $(FFLAGS) forward_solver.f90 

aniprop_subroutines.o : aniprop_subroutines.f
	$(FC) $(FFLAGS)  aniprop_subroutines.f

rf_aniso_subroutines.o : rf_aniso_subroutines.f
	$(FC) $(FFLAGS) rf_aniso_subroutines.f

util_subroutines.o : util_subroutines.f 
	$(FC) $(FFLAGS) util_subroutines.f 

eispack.o : eispack.f
	$(FC) $(FFLAGS) eispack.f

refft.o : refft.c
	$(CC) $(FFLAGS)  refft.c

clean :
	rm -f $(TARGETS) $(OBJS) 
