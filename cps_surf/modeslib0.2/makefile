include flags.mk
#for some reason, if the target (e.g. nms) as the same name
#as the directory (e.g. nms), gmake find that the target is
#"up to date" whatever it is the case or not. That is why I have
#added an "_" to the target (e.g. nms_)
all: minosy nms_ stock_boule_
minosy:
	cd yannos; $(MAKE) minosy
yannos_MPI: 
	cd yannos; $(MAKE) yannos_MPI
nms_:
	cd nms; $(MAKE) all
stock_boule_:
	cd stock_boule; $(MAKE) all
clean:
	cd yannos; $(MAKE) clean_minosy; $(MAKE) clean_yannos_MPI
	cd nms; $(MAKE) clean
	cd stock_boule; $(MAKE) clean
very_clean:
	$(RM) bin/*
	make clean
