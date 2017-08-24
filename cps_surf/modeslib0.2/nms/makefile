include ../flags.mk
LOBJ    =obj/
#
DBGFLAG = NO
FFLAGS = $(Opt) -Imod  $(MODOUT) 
LDFLAGS =  $(Opt) 
MYBIN=../bin

LDFLAGS = $(Opt) 

SRCS = nms.f90 

OBJS =	$(LOBJ)xerbla.o $(LOBJ)def_gparam.o  $(LOBJ)lsame.o $(LOBJ)zgemv.o \
        $(LOBJ)dfour1.o $(LOBJ)nrtype.o $(LOBJ)nrutil.o $(LOBJ)global_main.o $(LOBJ)time_function.o \
        $(LOBJ)module_modele.o $(LOBJ)module_spline.o   \
        $(LOBJ)util_fctp.o $(LOBJ)module_modes.o $(LOBJ)nms.o 

OBJS2 = $(LOBJ)def_gparam.o $(LOBJ)nrtype.o $(LOBJ)nrutil.o  $(LOBJ)global_main.o  \
        $(LOBJ)module_modele.o $(LOBJ)module_spline.o  $(LOBJ)util_fctp.o $(LOBJ)get_fctp.o

MODS =	

all: nms get_fctp generate_prem

nms: $(MYBIN)/nms
$(MYBIN)/nms: $(OBJS)
	$(F90)   $(LDFLAGS) $(OBJS) -o $(MYBIN)/nms 
	$(MOVE_MOD)
clean:
	\rm $(LOBJ)/*   mod/*

get_fctp: $(MYBIN)/get_fctp
$(MYBIN)/get_fctp: $(OBJS2)
	$(F90)   $(LDFLAGS) $(OBJS2) -o $(MYBIN)/get_fctp
	$(MOVE_MOD)

generate_prem: $(MYBIN)/generate_prem
$(MYBIN)/generate_prem: generate_prem.f90
	$(F90)   $(LDFLAGS)  generate_prem.f90   -o $(MYBIN)/generate_prem

$(LOBJ)dfour1.o: dfour1.f
	$(F90)   $(Opt9077) -c dfour1.f
	$(MV)  dfour1.o $(LOBJ)
$(LOBJ)zgemv.o: zgemv.f
	$(F90)   $(Opt9077) -c zgemv.f
	$(MV)  zgemv.o $(LOBJ)
$(LOBJ)lsame.o: lsame.f
	$(F90)   $(Opt9077) -c lsame.f
	$(MV)  lsame.o $(LOBJ)
$(LOBJ)xerbla.o: xerbla.f
	$(F90)   $(Opt9077) -c xerbla.f
	$(MV)  xerbla.o $(LOBJ)
$(LOBJ)%.o: %.f90 	
	$(F90)   $(FFLAGS) -c $*.f90
	$(MV)  $*.o $(LOBJ)
