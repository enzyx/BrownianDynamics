# BrownianDynamics
Developed by Fabian Zeller and Manuel P. Luitz(2017).

## Introduction
BrownianDynamics (BD) is a software package to simulate association events of ligand molecules to their receptor with Brownian Dynamics. The receptor and ligand are thereby treated as rigid molecule. The receptor is kept fixed at the origin of the coordinate system while the ligand can diffuse around it. The diffusive motions of the ligand account for the relative translational and rotational diffusion of both, receptor and ligand. To account for the long-range electrostatic interactions between ligand and receptor, the electrostatic field around the receptor has to be calculated first with the [Adaptive Poisson Boltzmann Solver](http://www.poissonboltzmann.org/docs/apbs-faq/). At the beginning of a BD simulation, the ligand is randomly positioned on a shell around the receptor. The simulation is terminated when either the ligand diffuses into the binding region (fulfills the binding criterion) or when the ligand moves away from the receptor by more than a threshold distance. Final ligand/receptor configurations of rajectories that reach the binding criterion are stored. From the fraction of trajectories where the ligand reaches the binding site the binding rate ("on" rate) of ligand and receptor can be estimated. 
The binding criterium can be either defined as center of mass (COM) distances between selected ligand/receptor atoms
or as the distance root mean square deviation (drmsd) between selected ligand/receptor atoms. The latter allows not only to define the interactive patch on the receptor surface but also to define the orientation between ligand and receptor for a successful binding event.

## Installation
You need the tools `apbs` and `pdb2pqr` which can be installed on a Ubuntu system with:

    $ apt-get install apbs pdb2pqr

BD itself will need `parmed` and `numpy` which are available with python pip:

    $ pip install numpy parmed

## HowTo
Assume we have two `pdb` files `receptor.pdb` and `ligand.pdb` in our working directory and we want to calculate the association rate of the ligand into a roughly define reaction surface around the receptors binding site.
Before we can perform the BD run we need to prepare `pqr` files from the `pdb` files that include partial charges for each atom. We can do this with the `pdb2pqr` tool:

$ pdb2pqr --ff amber receptor.pdb receptor.pqr
$ pdb2pqr --ff amber ligand.pdb ligand.pqr

The `--ff amber` option specifies that we want to use partial charges from the amber forcefield.
Next we need to calculate the electrostatic field by solving the [Poisson Boltzmann equation](https://en.wikipedia.org/wiki/Poisson%E2%80%93Boltzmann_equation) around the receptor.
We use the APBS software package to perform the PB calculations. A sample `apbs.in` configuration file is provide with the BD sources, the documentation can be found at [here](http://www.poissonboltzmann.org/docs/apbs-faq/). The sample config should normally work for you just recheck the input file path for the `receptor.pqr` so that `apbs` can find it. Run the `apbs` program as follows:

$ apbs apbs.in

This will generate a file `pot.dx` containing the electrostatic potential on a 3D grid.
After having calculated the electrostatic field, we can start the BD simulation.

$ ./BrownianDynamics/BrownianDynamics.py -c template.conf

Finally we can calculate the "on" rate from bulk state to the reaction surface ([S. H. Northrup et al., J. Chem. Phys. 80 (4), 15 February 1984](http://dx.doi.org/10.1063/1.446900)).

$ ./BrownianDynamics/BD_rate.py -c bd.conf


## Configuration file
A sample configuration file `template.conf` is provided as a starting point. The syntax is based on the python module `ConfigParser`. 
