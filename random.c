#include "stdlib.h"
#include "stdio.h"
#include "omp.h"

void main() {
 	 struct drand48_data drand_buf;
	 int seed;
	 double x,y;


#pragma omp parallel private(x, seed, drand_buf)
	 {
	int i;
	seed = 1202107158 + omp_get_thread_num() * 1999;
	srand48_r (seed, &drand_buf);

for (i=0; i<3; i++) {
	 
	 drand48_r (&drand_buf, &x);

	 printf("thread %d : %.3f\n",omp_get_thread_num(),x);
	
	#pragma omp barrier  

#pragma omp single

	 printf("\n");
}
	 }
}
