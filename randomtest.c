#include "stdio.h"
#include "ziggurat_openmp.c"

void main() {
  float fn[128];
  uint32_t kn[128];
  float wn[128];
	r4_nor_setup ( kn, fn, wn );
  	
  	uint32_t jsr=123456;

	int i;
	
	uint32_t seed[2];

	for (i=0; i<2; i++) 
		seed[i]=shr3_seeded(&jsr);

//	for (i=0; i<10; i++)  {
//		printf("%d\t%d\t", seed[0], seed[1]);
//		printf("%d\t%d\n", shr3_seeded(&seed[0]), shr3_seeded(&seed[1]));
//	}
//

	for (i=0;i<1000000; i++)
		printf("%f\n", r4_nor(&jsr,kn,fn,wn));

	

}
