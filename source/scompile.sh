#!/bin/bash

base="/global/software/sl-7.x86_64/modules/gcc/"
version="5.4.0/"

netcdfinclude=$base$version"netcdf/4.4.1.1-gcc-p/include"
netcdflib=$base$version"netcdf/4.4.1.1-gcc-p/lib"

fftwinclude=$base$version"fftw/3.3.6-gcc/include/"
fftwlib=$base$version"fftw/3.3.6-gcc/lib"

mpiinclude=$base$version"openmpi/2.0.2-gcc/include"
mpilib=$base$version"openmpi/2.0.2-gcc/lib"

#hdf5include=$base$version"hdf5/1.8.18-gcc-p/include"
#hdf5lib=$base$version"hdf5/1.8.18-gcc-p/lib"

localinclude="../include"

cd /global/home/users/michaeldiaz/Model/source

mpic++ -O3 -g -o ../solve.exe main.cpp mpidriver.cpp ensemble.cpp pcomm.cpp solver.cpp temperature.cpp initializer.cpp data_initializer.cpp interpolate.cpp fluxes.cpp surface.cpp pressure.cpp advection.cpp damping.cpp boundaries.cpp files.cpp util.cpp Heating.cpp microphysics.cpp kessler.cpp rutledge.cpp energy.cpp trajectory.cpp laplacian.cpp -L $fftwlib -lfftw3 -L $netcdflib -lnetcdf -L $mpilib -lmpi -lm -I $fftwinclude -I $netcdfinclude -I $mpiinclude -I $localinclude