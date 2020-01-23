#!/bin/bash

#------------------------------------------------
# Change to user specific paths
#------------------------------------------------
netcdf_path="/usr/local"
fftw_path="/usr/local"
mpi_path="/usr/local"

#------------------------------------------------
# Probably don't need to change below this line
#------------------------------------------------
netcdfinclude=$netcdf_path"/include"
netcdflib=$netcdf_path"/lib"

fftwinclude=$fftw_path"/include"
fftwlib=$fftw_path"/lib"

mpiinclude=$mpi_path"/include"
mpilib=$mpi_path"/lib"

localinclude="../include"

cd source

mpic++ -O3 -g -o ../bin/solve.exe main.cpp mpidriver.cpp solver.cpp ensemble.cpp initializer.cpp data_initializer.cpp interpolate.cpp fluxes.cpp pressure.cpp advection.cpp surface.cpp damping.cpp boundaries.cpp pcomm.cpp files.cpp util.cpp Heating.cpp kessler.cpp rutledge.cpp microphysics.cpp energy.cpp temperature.cpp trajectory.cpp laplacian.cpp process_input.cpp -L $fftwlib -lfftw3 -L $netcdflib -lnetcdf -L $mpilib -lmpi -lm -I $fftwinclude -I $netcdfinclude -I $mpiinclude -I $localinclude

#mpicc -O3 -g -o ../bin/solve.exe main.cpp mpidriver.cpp solver.cpp ensemble.cpp initializer.cpp data_initializer.cpp interpolate.cpp fluxes.cpp pressure.cpp advection.cpp surface.cpp damping.cpp boundaries.cpp pcomm.cpp files.cpp util.cpp Heating.cpp kessler.cpp rutledge.cpp microphysics.cpp energy.cpp temperature.cpp trajectory.cpp laplacian.cpp process_input.cpp -L $libpath -lnetcdf -lmpi -lfftw3 -lm -I $includepath -I $localinclude