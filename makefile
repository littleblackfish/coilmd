somon: avx-single

avx-single : src/dna.c
	icc -march=core-avx-i -Ofast -DCIRCULAR dna.c -o cir-single
	icc -march=core-avx-i -Ofast  dna.c -o lin-single

gcc-single : src/dna.c
	gcc -lm -O3 src/dna.c -o lin
	gcc -lm -O3 -DCIRCULAR src/dna.c -o cir

gcc-multi : src/dna.c
	gcc -fopenmp -lm -O3 src/dna.c -o lin
	gcc -fopenmp -lm -O3 -DCIRCULAR src/dna.c -o cir


	
