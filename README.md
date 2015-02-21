# coilmd
MD simulation for a coarse grained DNA model, parallelized by openMP

## Compilation

icc -Ofast -openmp dna.c -o coil 

gcc -O3 -fopenmp -lm dna.c -o coil

It also compiles fine without OpenMP support. ICC will throw some warnings but they are fine. 

You can also define NDEBUG (-DNDEBUG) to output neighbour counts to a seperate file.

## Coil, ladder, circle

The program compiles to simulate a linear coil by default. 

you can add *-DLADDER* to simulate a flat ladder structure (without any twist associated with it) or *-DCIRCULAR* to simulate a circualr coil.

## Running

The final executable will take 3 parameters, N, nsteps, temperature, in this order. 

N is the base count, not particle count, so there will be 2N particles in total. 

If it finds a compatible restart file in the working directory, it will pick up from there. 

If it finds an incompatible restart file it will halt to avoid possibly overwriting something. 

## Output

The program will normally output 3 files, energy.dat, traj.vtf. bubbles.dat, which are the energy and position trajectories of the system and a huge matrix that keeps the binding state of each 'base pair'. 

### traj.vtf

This is a MD trajectory file directly readable by VMD. It is much less efficient than for example dcd but it does the job. A great deal of space can be saved by converting it to dcd using vmd, but then a seperate topology file is necessary and the dcd cannot be read easily as it is binary. The first 2N lines of the vtf file can be used as a vsf file to keep the topology only. 

The vmd scripts to get this done are work in progress, to be included in the repo eventually. 

### energy.dat

This is a file that keeps some energy sums for each timestep for which the coordinates are printed for. The exact order is :

* timestep
* temperature
* intra-strand bond energy
* inter-strand bond energy
* dihedral energy
* hardcore repulsive (non-bonded) energy

### bubbles.dat

This is a binary matrix of size *nsteps x N*. Each line represents a snapshot, each 0 represents a broken inter-strand bond and each 1 represents an intact inter-strand bond. Majority of our calculations are based on this information only. 

Note that this file may grow very big but also has tremendous potential for compression. 
