#!/bin/bash

#/opt/cray/pe/fftw/3.3.8.2/haswell/

netcdfinclude="/opt/cray/pe/netcdf/4.6.3.0/intel/19.0/include"
netcdflib="/opt/cray/pe/netcdf/4.6.3.0/intel/19.0/lib"

#netcdfinclude="/opt/cray/pe/netcdf-hdf5parallel/4.6.3.0/intel/19.0/include"
#netcdflib="/opt/cray/pe/netcdf-hdf5parallel/4.6.3.0/intel/19.0/lib"

fftwinclude="/opt/cray/pe/fftw/3.3.8.2/haswell/include/"
fftwlib="/opt/cray/pe/fftw/3.3.8.2/haswell/lib"

mpiinclude="/global/common/software/m3169/cori/openmpi/3.1.4/intel/include"
mpilib="/global/common/software/m3169/cori/openmpi/3.1.4/intel/lib"

#hdf5include="/opt/cray/pe/hdf5-parallel/default/intel/16.0/include"
#hdf5lib="/opt/cray/pe/hdf5-parallel/default/intel/16.0/lib"

localinclude="../include"
#export NETCDF_DIR=/opt/cray/pe/netcdf-hdf5parallel/default/intel/16.0
#export PATH=/opt/cray/pe/netcdf-hdf5parallel/default/bin:$PATH

cd /global/homes/m/michaeld/Model/source

mpic++ -O3 -g -o ../solve.exe main.cpp mpidriver.cpp ensemble.cpp pcomm.cpp solver.cpp temperature.cpp initializer.cpp data_initializer.cpp interpolate.cpp fluxes.cpp surface.cpp pressure.cpp advection.cpp damping.cpp boundaries.cpp files.cpp util.cpp Heating.cpp microphysics.cpp kessler.cpp rutledge.cpp energy.cpp trajectory.cpp laplacian.cpp -L $fftwlib -lfftw3 -L $netcdflib -lnetcdf -L $mpilib -lmpi -lm -I $fftwinclude -I $netcdfinclude -I $mpiinclude -I $localinclude