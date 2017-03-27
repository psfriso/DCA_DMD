include ../config.mk
#
OBJS = IO.o random.o energy.o colisio.o DMD.o
MODS = paramSet.o  Structure.o ls_rmsd.o sort.o stepPotentials.o evolution.o simulation.o
EXEFILE = discrete

discrete: ${MODS} $(OBJS) $(LIBDIR) main.o
	$(F77) $(LFLAGS) -fopenmp $(ENDFLAG) -o $(EXEFILE) $(LOPTIONS) $(OBJS) $(MODS) main.o -l Discrete 

main.o: vars.h

Structure.o: Structure.f 
	$(F77) -ffree-form $(OPTIONS) $(CFLAGSNOCHECK) Structure.f

.f.o: $<
	$(F77) -ffree-form $(OPTIONS) $(CFLAGS) $<
.F.o: $<
	$(F77) -ffree-form $(OPTIONS) $(CFLAGS) $<
clean: 
	rm *.o *.mod $(EXEFILE)
