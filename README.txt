This program is an idealized numerical model that solves the anelastic equations of motion for the atmosphere based on a perturbation / basic state decomposition. It includes a microphysics and turbulence parameterization and can be run on parallel architecture. The equations and the numerical methods used to solve them are described in more detail in Diaz and Boos (2019a) and Diaz and Boos (2019b).

-------------------------------------------------------------
I. COMPILATION
-------------------------------------------------------------

This program requires three external libraries: NetCDF (the file format, downloadable at https://www.unidata.ucar.edu/software/netcdf/), FFTW3 (a Fast Fourier Transform library, downloadable at http://www.fftw.org), and an implementation of MPI (the parallel programming library, for example, OpenMPI). To compile, locate the file "compile.sh". In this file, you will need to set the file paths to your own compiled versions of the NetCDF, FFTW3, and MPI libraries. The directories pointed to by these paths must contain the "include" and "lib" subdirectories to each library. The MPI compiler wrappers (e.g. mpic++) must also be in the execution path. Run the script with "./compile.sh". If compilation is successful, the executable "solve.exe" will appear in the "bin" folder.

To switch between the parallel and serial version, locate the file "stdafx.h" in the "include" directory. Locate the line with "#define PARALLEL". Set it to "#define PARALLEL 0" for the serial version and "#define PARALLEL 1" for the parallel version. The program must be recompiled when switching between parallel and serial.

-------------------------------------------------------------
II. RUNNING
-------------------------------------------------------------

--------------------------------------------
A. Command-Line Options and Model Settings
--------------------------------------------

Before running the model, you will probably need to provide input data (see section B below).

The program is run using a specially formatted input file, which controls the model settings (grid spacing, grid size, model physics, etc.). The original version of the file is named "input_params.txt". Open it and change the values accordingly. To run, type the following line into the terminal window:

mpirun -np 4 ./solve.exe -f input_params.txt 

The number following "-np" is the number of processes to use. If this doesn't work, you may need to set the execution path to include the "bin" file in the MPI directory.

A basic state with both horizontal and vertical shear from a soon to be submitted paper can be generated using the command line option "-s". For example,

mpirun -np 4 ./solve.exe -f input_params.txt -s "1.0 -2.0 -1.0"

where the first value is the maximum vorticity (x 10^-4 1/s), the second is the upper-level shear (x 10^-3 m/s/km), and the third is the lower-level shear (x 10^-3 m/s/km). In order for this option to work, the value for "basic_state_init_option" needs to be set to "1" to direct the program to generate the basic state from code and not an input file.

--------------------------------------------
B. Model Input
--------------------------------------------

There are two ways to generate model input: either from an external NetCDF file or from the program itself. These options are set after "basic_state_input_file" and "perturbation_input_file" in the input text file, with "0" being from an external file and "1" being from code. Before running, you'll need to ensure that the model either has an input NetCDF file or will generate its initial perturbation and/or basic state from code. The perturbation and basic state can be set to the same file name if the basic state and perturbation are in the same file. The input grid size and spacing can be the same as or different from that at which the model is run.

Currently, the only option to generate a basic state from code is the one listed in the previous section (i.e. command-line option "-s"). The only option to generate an initial perturbation from code is an axisymmetric vortex based on Murthy and Boos (2018). To use this option, find the section labeled "&initialization_constants", set the values as desired, and make sure that "perturbation_init_option" is set to "1" to initialize from code and not an input NetCDF file.

The python scripts "sample_initializations.py" and "initialize.py" are useful for generating input files. The function initialize() converts arrays of meteorological data into a NetCDF file formatted to be read as input into the model. The file "sample_initializations.py" has a number of sample initializations that demonstrate how to use this function.

--------------------------------------------
C. Model Output
--------------------------------------------

The model output is given as "output_file" in the input text file. The variables in the output file are as follows:

u-wind(t, x, y, z)	-> perturbation zonal wind (m/s)
v-wind(t, x, y, z)	-> perturbation meridional wind (m/s)
w-wind(t, x, y, z)	-> perturbation vertical wind (m/s) 
theta(t, x, y, z) 	-> perturbation potential temperature (K)
pi(t, x, y, z)	 	-> perturbation pressure (Pa / density)
int_fric(t, x, y) 	-> vertically integrated friction
rainfall(t, x, y)	-> accumulated rain fall (mm)
snowfall(t, x, y)	-> accumulated snow fall (mm)
qv(t, x, y, z)  	-> perturbation water vapor mixing ratio (g/g)
qc(t, x, y, z)   	-> cloud water mixing ratio (g/g)
qr(t, x, y, z)   	-> rain water mixing ratio (g/g)
qi(t, x, y, z)   	-> ice water vapor mixing ratio (g/g)
qs(t, x, y, z) 		-> snow water vapor mixing ratio (g/g)
fric(t, x, y, z) 	-> loss of kinetic energy from turbulence parameterization  (J)
ubar(x, y, z) 		-> basic state zonal wind (m/s)
vbar(x, y, z) 		-> basic state meridional wind (m/s)
wbar(x, y, z) 		-> basic state vertical wind (m/s)
thbar(x, y, z)		-> basic state potential temperature (K)
qbar(x, y, z)		-> basic state water vapor mixing ratio (g/g)
pbar(x, y, z)		-> basic state pressure (non-dimensional through Exner function)
tb(z)				-> base state potential temperature (K)
pib(z) 				-> base state pressure (non-dimensional through Exner function)
qb(z) 				-> base state water vapor mixing ratio (g/g)
zu(z)				-> height of scalar variables (m)

Note that most of the variables are split into perturbation, basic state, and base state and will need to be added together to get the full state. This could require adding together two (for winds) or three (for scalars) fields. One exception is the condensate variables, namely mixing ratios of rain (qr), snow (qs), ice (qi), and cloud (qc), which don't have a basic state. In contrast with "tbar" and "qbar", the pressure variable "pbar" represents the full basic state, with "pib" already included. The pressure variables are also different in that the perturbation and basic states have different units.

The model output files are formatted similarly to the input files and can therefore be read as input. The input text file provides several ways to exploit this similarity using the options "is_restart_run", "perturbationFileTime", "startOutfileAt", and "create_new_output_file". However, be aware that restarting from a model output file may yield slightly different results than running straight through, since the output data is stored as single precision floating point values but the model is run with double precision.

The model can also output multiple budgets, including temperature, moisture, potential energy, potential vorticity, and stream function. Note that these budgets are accumulated between output times and are therefore not instantaneous. Currently, the budget work for only the parallel version of the code.

-------------------------------------------------------------
III. PLOTTING
-------------------------------------------------------------

A sample python script named "plot_model_output.py" is provided to plot the model output.
