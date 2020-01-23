
CC=mpicc
IDIR=/usr/local/include
LDIR=/usr/local/lib


LOCALIDIR=../include
CFLAGS=-O3 -I$(IDIR) -I$(LOCALIDIR)



ODIR=obj
LIBS = -lnetcdf -lmpi -lfftw3 -lm

_OBJ = main.o mpidriver.o solver.o ensemble.o initializer.o data_initializer.o interpolate.o fluxes.o pressure.o advection.o surface.o damping.o boundaries.o pcomm.o files.o util.o Heating.o kessler.o rutledge.o microphysics.o energy.o temperature.o trajectory.o laplacian.o process_input.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


solvemake: $(OBJ)
	$(CC) $(CFLAGS) $(LIBS) -o ../bin/solve.exe ${OBJ}

_DEPS = stdafx.h energy.h initializer.h trajectory.h temperature.h ensemble.h data_initializer.h process_input.h
DEPS = $(patsubst %,$(LOCALIDIR)/%,$(_DEPS))

$(ODIR)/main.o: main.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

#-L $libpath -lnetcdf -lmpi -lfftw3 -lm -I $includepath -I $localinclude


#mpicc -O3 -g -o ../bin/solve.exe

clean:
	rm $(ODIR)/*.o