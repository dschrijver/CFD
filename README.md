# Incompressible Navier-Stokes solver
## Requirements and Installation
On Ubuntu, base requirements can be installed using the following command,

```
sudo apt install binutils openmpi-bin libopenmpi-dev zlib1g zlib1g-dev make g++
```

I use HDF5 for data output. If not yet present on your device, HDF5 can be installed in the root folder of this directory using the following command,

```
make install_hdf5
```
## Running simulation
Settings for the simulation can be set in [params.h](params.h). Note that the variable poisson_epsilon denotes the maximal allowed divergence in the velocity field. Initial conditions and boundary conditions can be chosen in [definitions.h](definitions.h) by uncommenting and commenting lines. Non-zero velocities can be imposed at boundaries by uncommenting these lines and changing the desired velocity following the macro. 

After setting up both [params.h](params.h) and [definitions.h](definitions.h), you can start compilation of the program. If you want to make use of an HDF5 build that was not installed using ```make install_hdf5``` into the root directory of this repository, you have to change the first two lines in [Makefile](Makefile). Any additional optimization flags, for instance the architecture-specific march flag, can be added to the compilation by uncommenting the third line in [Makefile](Makefile), and possibly adding additional flags.

The code can then be compiled with either ```make cfd``` or just ```make```. 
## Testing Poisson solver
The Poisson solver can be tested in the folder [test_poisson](test_poisson/). Settings can be changed again in [test_poisson/params.h](test_poisson/params.h) and [test_poisson/definitions.h](test_poisson/definitions.h). The flag INI_RHS_SINE needs to be combined with the periodic boundary conditions. The flag INI_RHS_COSINE needs to be combined with the no-slip boundary conditions. The flag INI_RHS_NOISE can be tested with any set of boundary conditions, and is used to stress-test the Poisson solver. 

The test can be compiled with ```make```. The INI_RHS_SINE can be validated, after running the test, using ```python ../scripts/validate_periodic.py```. INI_RHS_COSINE can be validated, after running the test, using ```python ../scripts/validate_noslip.py```. Finally, the results for INI_RHS_NOISE can be visualized using ```python ../scripts/plot_noise_solution.py```.