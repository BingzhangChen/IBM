.SUFFIXES: .o .i .f90 .f
.f90.o:
	$(CFT) -c $(FFLAGS) $*.f90 -o $*.o
.f.o:
	$(CFT) -c $(FFLAGS) $*.f -o $*.o

MAIN = IBM
#Source file
SRCS = Readcsv.f90 tridiagonal.f90\
	multiGauss.f90\
	Trait_functions.f90\
    Calc_PAR.f90\
	params.f90\
	grid.f90 gridinterp.f90 time_interp.f90\
	compute_Kv.f90\
	variables.f90\
	Time_settings.f90\
	forcing.f90\
	lagrange.f90\
	IO_files.f90\
 	Advection_center.f90 Diff_center.f90\
	Geider_Lag.f90\
	initialization.f90\
	timestep.f90\
	Main.f90

OBJS =$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))

#Executable file
$(MAIN): $(OBJS)
	$(CFT) $(FFLAGS) $(LFLAGS) -o $(MAIN) $(OBJS)
