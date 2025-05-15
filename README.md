# Gizmo Pipeline

This is a pipeline to create trace particle datasets from Gizmo simulation code. There are main steps to create the dataset:
1. MakeCloud: Creates an initial condition for the simulation and specifies the parameters for the Gizmo simulations.
2. Gizmo/STARFORGE: Configure and compile the Gizmo code with the desired physics and then run the simulation.
3. Trace The Cells: Trace the cells in the simulation to create the dataset for training.

## Code Used
- [MakeCloud](https://github.com/mikegrudic/MakeCloud/tree/master)
- [Gizmo](https://github.com/pfhopkins/gizmo-public)

## MakeCloud
MakeCloud creates the initial conditions and the parameter file for the simulations. It is written python and uses the following libraries:
- numpy
- h5py
- scipy
- docopt

Additionally, you will need to install the glassy initial conditions and place it in the MakeCloud directory. You can download the glassy initial conditions from [here](https://data.obs.carnegiescience.edu/starforge/glass_orig.npy).

In the MakeCloud directory, there are two shells scripts to setup the initial conditions for gravitational collapse and for turbulent MHD. 

## STARFORGE_FILES
The STARFORGE_FILES directory contains the files used to run the simulations with STARFORGE. Due to different versions of Intel installed on each machine, this can only run on Frontera using intel/19. To compile and run a STARFORGE simulation follow these steps:
1. Uncomment the "Frontera" system in the `Makefile.systype` file. 
1. Change the `Makefile` variable `CONFIG` to your desired configuration file. The directory `STARFORGE_FILES` contains the following configuration files:
    - `Grav_rad_Config.sh`: This is the configuration file for running gravitational collapse with radiation.
    - `Config_Rad.sh`: This is the configuration file for running a turbulent hydrodynamics with radiation and no gravity.
    - `Config_Rad_MHD.sh`: This is the configuration file for running a turbulent magneto-hydrodynamics with radiation and no gravity.
1. On a development node, load the proper modules required to run STARFORGE:
    ```
    module purge
    module load TACC intel impi hdf5 gsl fftw3
    ```
1. Compile the code:
    ```
    make clean
    make
    ```
1. Run the code:
    ```
    export OMP_NUM_THREADS=2
    ibrun ./GIZMO <paramfile> 1>output.txt 2>error.txt
    ```
    where `<paramfile>` is the parameter file created by MakeCloud. In STARFORGE_FILES, I included two example parameter files for the gravitational collapse, `params_M5.9_R0.025392_Z1_S0_A0_B0_I1_Res22_n2_sol0.5_42.txt` and turbulent/MHD,`params_M6e2_R1.5_Z1_S0_A2_B0.01_I1_Res32_n2_sol0.5_42.txt`. 
## Trace The Cells
The Trace The Cells code is used to track the chemistry input parameters from the simulations. It is written in python and uses the following libraries:
- numpy
- h5py
- glob
- matplotlib
- pytreegrav
The input parameters for the trace cell script are at the top of the `trace_cells.py` file. The input parameters are:
    - `bands`: The number of radiation bands in the simulation.
    - `path`: The directory where the simulation is located.
    - `radius`: The radius of the sphere in the simulation. All particles within this radius will be traced.
    - `rays`: The number of rays used to calculate the number density of the particles.
    - `ray_saved`: Which of the rays to save in the dataset.

Running the script generates a dataset with the parameters described in the dataset readme.