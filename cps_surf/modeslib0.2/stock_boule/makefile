include ../flags.mk
MOD=mod/
LOBJ=obj/
BIN = ../bin

FFLAGS = $(Opt) -I$(MOD) $(MODOUT)  
FFLAGS7790 = $(Opt9077) -I$(MOD) $(MODOUT) $(INC)
LDFLAGS = $(Opt) 

obj= $(LOBJ)beschb.o  $(LOBJ)bessjy.o $(LOBJ)sphbes.o \
         $(LOBJ)chebev.o $(LOBJ)bsf.o $(LOBJ)gauleg.o     \
         $(LOBJ)def_gparam.o $(LOBJ)module_boule.o $(LOBJ)stock_boule.o
all: stock_boule
stock_boule:  $(BIN)/stock_boule
$(BIN)/stock_boule: $(obj)
	$(F90)  $(FFLAGS) $(obj) -o $(BIN)/stock_boule

$(LOBJ)beschb.o: beschb.f
	$(F90)  $(FFLAGS7790)   -c beschb.f 
	$(MV)  beschb.o $(LOBJ)
$(LOBJ)bessjy.o: bessjy.f
	$(F90)  $(FFLAGS7790)   -c bessjy.f 
	$(MV)  bessjy.o $(LOBJ)
$(LOBJ)sphbes.o: sphbes.f
	$(F90)  $(FFLAGS7790)   -c sphbes.f 
	$(MV)  sphbes.o $(LOBJ)
$(LOBJ)chebev.o: chebev.f
	$(F90)  $(FFLAGS7790)   -c chebev.f 
	$(MV)  chebev.o $(LOBJ)
$(LOBJ)bsf.o: bsf.f
	$(F90)  $(FFLAGS7790)   -c bsf.f 
	$(MV)  bsf.o $(LOBJ)
$(LOBJ)gauleg.o: gauleg.f
	$(F90)  $(FFLAGS7790)   -c gauleg.f 
	$(MV)  gauleg.o $(LOBJ)
$(LOBJ)%.o: %.f90 	
	$(F90)   $(FFLAGS) -c $*.f90
	$(MV)  $*.o $(LOBJ)

clean: 
	\rm $(obj) $(MOD)/*