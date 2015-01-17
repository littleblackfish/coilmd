# coilmd
MD simulation for a coarse grained DNA model, parallelized by openMP

compile with 

icc -Ofast -openmp dna.c -o dna
gcc -O3 -fopenmp -lm dna.c -o dna
