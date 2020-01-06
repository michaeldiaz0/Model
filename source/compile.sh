#!/bin/bash

includepath="/usr/local/include"
localinclude="../include"
libpath="/usr/local/lib"

cd /Users/michaeldiaz/Desktop/Model/source

mpicc -O3 -g -o ../bin/solve.exe main.cpp mpidriver.cpp solver.cpp ensemble.cpp initializer.cpp data_initializer.cpp interpolate.cpp fluxes.cpp pressure.cpp advection.cpp surface.cpp damping.cpp boundaries.cpp pcomm.cpp files.cpp util.cpp Heating.cpp kessler.cpp rutledge.cpp microphysics.cpp energy.cpp temperature.cpp trajectory.cpp laplacian.cpp process_input.cpp -L $libpath -lnetcdf -lmpi -lfftw3 -lm -I $includepath -I $localinclude
