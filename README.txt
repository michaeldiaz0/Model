This is idealized model which solves the equations of motion for the atmosphere for a perturbation / basic state decomposition. It is described in more detail in Diaz and Boos (2019a) and Diaz and Boos (2019b).

I. COMPILATION

To compile, locate the file "compile.sh". In this file, you will need to set the file paths to your own NETCDF, FFTW3, and MPI libraries. These paths must contain the "include" and "lib" directories to each library. Run the script with "./compile.sh".

To switch between the parallel and serial version, locate the file "stdafx.h" in the "include" directory. Locate the line with "#define PARALLEL". Set it to "#define PARALLEL 0" for the serial version and "#define PARALLEL 1" for the parallel version. The program must be recompiled when switching between parallel and serial.

II. RUNNING

A. Command line options and model settings

The program is run using a specially formatted input file, which controls the model settings. The original version of the file is called "input_params.txt". Open it and change the values accordingly. To run, type the following line into the terminal window:

mpirun -np 4 ./solve.exe -f input_params.txt 

The number following "-np" is the number of processes to use.

The basic state from paper 3 can be generated using the command line option "-s". For example,

mpirun -np 4 ./solve.exe -f input_params.txt -s "1.0 -2.0 -1.0"

where the first value is the vorticity (x 10^-4 1/s), the second is the upper-level shear (x 10^-3 m/s/km), and the third is the lower-level shear (x 10^-3 m/s/km). In order for this option to work, the value for "basic_state_init_option" needs to be set to "1" to direct the program to generate the basic state from code and not an input file.

B. Model input and output files

Before running, you'll need to ensure that the model has an input NETCDF file to provide the initial conditions and basic state. This is set on the lines containing "basic_state_input_file" and "perturbation_input_file" in the input text file. The former provides the basic state and the latter the perturbation. These can be set to the same file name if the basic state and perturbation are in the same file. The input grid size and spacing can be the same as or different from that at which the model is run. The model output is given as "output_file".

The python script "initialize.py" is useful for generating input files. The function initialize() converts arrays of meteorological data into a netcdf file formatted to be read as input into the model.

The model can generate an axisymmetric vortex based on Murthy and Boos (2018) for use as an initial perturbation. Find the section labeled "&initialization_constants", set the values as desired, and make sure that "perturbation_init_option" is set to "1" to initialize from code and not an input file.

III. PLOTTING

A sample python script is provided to plot the model output.
