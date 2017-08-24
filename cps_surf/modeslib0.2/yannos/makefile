include ../flags.mk
MOD=mod/
LOBJ=obj/
BIN = ../bin
#
FFLAGS = $(Opt) -I$(MOD) $(MODOUT)  
FFLAGS7790 = $(Opt9077) -I$(MOD) $(MODOUT) $(INC)
LDFLAGS = $(Opt) 

SRCS =	def_gparam.f90 yannos.f90 module_com_minos.f90  module_neominos.f90 \
	module_modele.f90 module_spline.f90  sub_minosa4.f util_minos.f90  \
	nrutil.f90  nrtype.f90 yannos_flag.f90 yannos_MPI.f90 layer.f90

OBJS =	$(LOBJ)def_gparam.o   $(LOBJ)module_com_minos.o $(LOBJ)module_modele.o $(LOBJ)util_minos.o \
        $(LOBJ)module_neominos.o  $(LOBJ)module_spline.o $(LOBJ)sub_minosa4.o\
        $(LOBJ)nrutil.o  $(LOBJ)nrtype.o   $(LOBJ)yannos_flag.o $(LOBJ)layer.o 

MODS =	$(MOD)module_com_minos.mod $(MOD)yannos_flag.mod\
        $(MOD)module_modele.mod $(MOD)module_neominos.mod \
        $(MOD)module_spline.mod $(MOD)util_minos.mod \
        $(MOD)nrutil.mod  $(MOD)nrtype.mod $(MOD)def_gparam.mod $(MOD)layer.mod 
#
all: minosy
minosy: $(BIN)/minosy
$(BIN)/minosy: $(OBJS)  $(LOBJ)yannos.o
	$(F90) $(LDFLAGS) $(LOBJ)yannos.o $(OBJS) $(LLIBS) -o $(BIN)/minosy 
	$(MOVE_MOD)
clean_minosy:
	\rm -f $(OBJS) mod/*  $(LOBJ)yannos.o   $(LOBJ)/*
#
yannos_MPI: $(BIN)/yannos_MPI
$(BIN)/yannos_MPI: $(OBJS) $(LOBJ)yannos_MPI.o
	$(MF90) $(LDFLAGS) -o $(BIN)/yannos_MPI  $(LOBJ)yannos_MPI.o $(OBJS) $(LLIBS)
	$(MOVE_MOD)
clean_yannos_MPI:
	\rm -f $(OBJS) $(LOBJ)yannos_MPI.o
clean:
	\rm $(OBJS) $(MODS)/*$(pM)   $(LOBJ)yannos.o   $(LOBJ)yannos_MPI.o
#
.SUFFIXES: $(SUFFIXES) .f90
$(LOBJ)yannos.o: $(LOBJ)yannos_flag.o $(LOBJ)module_minos.o $(LOBJ)module_neominos.o $(LOBJ)util_minos.o $(LOBJ)layer.o
$(LOBJ)yannos_MPI.o: $(LOBJ)yannos_flag.o $(LOBJ)module_minos.o $(LOBJ)module_neominos.o $(LOBJ)util_minos.o $(LOBJ)layer.o
$(LOBJ)module_minos.o: $(LOBJ)module_com_minos.o
$(LOBJ)sub_minosa4.o: $(LOBJ)module_com_minos.o $(LOBJ)yannos_flag.o
$(LOBJ)module_neominos.o: $(LOBJ)module_com_minos.o $(LOBJ)module_spline.o \
    $(LOBJ)module_modele.o $(LOBJ)util_minos.o $(LOBJ)yannos_flag.o
$(LOBJ)util_minos.o: $(LOBJ)yannos_flag.o $(LOBJ)layer.o
$(LOBJ)module_spline.o: $(LOBJ)nrtype.o $(LOBJ)nrutil.o   
#
$(LOBJ)sub_minosa4.o: sub_minosa4.f
	$(F90)   $(FFLAGS7790) -c sub_minosa4.f
	$(MV)  sub_minosa4.o $(LOBJ)
$(LOBJ)yannos_MPI.o: yannos_MPI.f90 	
	$(MF90)   $(FFLAGS) -c yannos_MPI.f90
	$(MV)  yannos_MPI.o $(LOBJ)
$(LOBJ)%.o: %.f90 	
	$(F90)   $(FFLAGS) -c $*.f90
	$(MV)  $*.o $(LOBJ)


