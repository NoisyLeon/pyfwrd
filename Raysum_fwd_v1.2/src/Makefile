# Makefile for raysum programs.
# Using f77 on Sun -- use -C -g when debugging, and -fast otherwise
# -autopar -parallel makes parallel code
FC = gfortran
#FFLAGS = -C -g -ftrap=%all,no%inexact
#LDFLAGS = -C -g
#FFLAGS = -fast
#LDFLAGS = -fast
# Using g77 on Sun or FreeBSD: (add -g for debugging)
# use -O -funroll-loops -m486
#F77 = g77
#FFLAGS = -Wall -O -funroll-loops -m486
#LDFLAGS =
# Using f2c on FreeBSD
#F77 = f77
#FFLAGS = -O -funroll-loops
#LDFLAGS =

EIG = eigenvec.o matrixops.o eispack-cg.o
RAYSUM = raysum.o $(EIG)
IO = readwrite.o
TRACE = phaselist.o buildmodel.o trace.o
MISFIT = misfit.o

# implicit rule for Fortran compilation. Changing the parameter header
# file recompiles everything.
.f.o: params.h
	$(FC) $(FFLAGS) -c $<
        
default: all

all: seis-spread seis-misfit

seis-spread: seis-spread.o $(RAYSUM) $(IO) $(TRACE)
	$(FC) $(LDFLAGS) -o seis-spread seis-spread.o $(RAYSUM) \
                             $(IO) $(TRACE)
	cp seis-spread ../bin
        
seis-misfit: seis-misfit.o $(MISFIT) $(IO)
	$(FC) $(LDFLAGS) -o seis-misfit seis-misfit.o $(MISFIT) $(IO)
	cp seis-misfit ../bin
clean:
	/bin/rm -f *.o *.core seis-spread seis-misfit

