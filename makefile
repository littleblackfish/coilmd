somon: avx-single

avx-single : avx-cir-single avx-lin-single

avx-cir-single : dna.c
	icc -march=core-avx-i -Ofast -DCIRCULAR dna.c -o cir-single

avx-lin-single : dna.c
	icc -march=core-avx-i -Ofast  dna.c -o lin-single
	
