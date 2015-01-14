# coilmd
MD simulation for a coarse grained DNA model, parallelized by openMP

compile with 

icc -Ofast -openmp -std=gnu99 dna.c -o dna
gcc -O3 -fopenmp -lm -std=gnu99 dna.c -o dna
