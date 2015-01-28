#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#ifdef _OPENMP
  #include "omp.h"
#else 
   int omp_get_thread_num() {return 0;}
   int omp_get_num_threads() {return 1;}
   int omp_get_max_threads() {return 1;}
#endif

#include "ziggurat_openmp.c"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//#define N 100
//#define NSTEPS 1000000
#define WFREQ 1000
#define DT 0.01
#define GAMMA 1

#define INTRA_BOND_LENGTH 0.5
#define INTER_BOND_LENGTH 1.1

#define INTER_BOND_CUT 1.2
#define HARD_CUT 0.4

#define K_BOND 100
#define MASS K_BOND

#define PHI_1 -100
#define PHI_2 -95

#define K_DIHEDRAL 1
#define EPSILON_DIHEDRAL 0.1

#define NEIGH_CUT 1.2

// (NEIGH_CUT/HARD_CUT)^3
#define MAX_NEIGH 64


static void printMat(float [][3]) ;
static void printBubble (FILE *) ;
static void zero(float [][3]) ;
static float calcTemp();
static float ziggurat(int num_thread) ;  
static float maxForce();

// global variables

float intraE, interE, dihedralE, hardE;

float x[2*N][3];
float xRef[2*N][3];
float v[2*N][3];
float f[2*N][3];
int neigh[2*N][MAX_NEIGH+1];
int isBound[N];

// used for ziggurat
static float fn[128];
static uint32_t kn[128];
static float wn[128];
// used for SHR3
static uint32_t *seed;

// some constants for performance 

static const float neighCutSq  = NEIGH_CUT * NEIGH_CUT ;
static const float neighSkinSq = (NEIGH_CUT-HARD_CUT) * (NEIGH_CUT - HARD_CUT) ;
static const float hardCutSq = HARD_CUT*HARD_CUT;

#include "restart.c"
#include "generators.c"
#include "neighbour.c"
#include "vtf.c"
#include "dihedral.c"
#include "harmonic.c"
#include "hardcore.c"
#include "harcos.c"
#include "langevin.c"

void main(int argc, char ** argv ) {

	#ifdef _OPENMP
		printf("Compiled with OpenMP, running %d threads.\n", omp_get_max_threads());
	#else
		printf("Compiled without OpenMP.\n");
	#endif
	
	if (argc<3)  {
		printf("I cannot run without nsteps.\n");
		exit(1);
	}
	if (argc<2)  {
		printf("I cannot run without temperature.\n");
		exit(1);
	}

	int nsteps = atoi(argv[1]);
	float temperature = atof (argv[2]);
	
	int t, i, rebuildCount = 0; 

	// seeding a random stream for each thread
	r4_nor_setup ( kn, fn, wn );
	uint32_t  jsr = 123456;
  	
	int max_threads = omp_get_max_threads();
	seed = ( uint32_t * ) malloc ( max_threads * sizeof ( uint32_t ) );
	
	for (i=0; i < max_threads; i++ ) 
		seed[i] = shr3_seeded ( &jsr );


	isBound[0]=1;
	isBound[N-1]=1;

	FILE *minim = initVTF("minim.vtf"); 
	FILE *traj  = initVTF("traj.vtf");
	FILE *energy  =	fopen("energy.dat", "w"); 
	FILE *neighCount   = fopen("neigh.dat", "w");
	FILE *bubbles =	fopen("bubbles.dat", "w");
	

	// minimization via Langevin at 0 temperature
	if ( !readRestart("restart") ) {
		genDNA(10.5);
		zero(f);
		zero(v);
		zero(xRef);

		printf ("Minimizing...");
		t=0;

		do {
			integrateLangevin(0.01,0);
			if (t%1000 ==0) {
				writeVTF(minim);
			}
			t++;
		} while (maxForce() > 0.1 && t <1000000 );

		printf ("done after %d steps.\n",t);
	}

	printf("N=%d, T=%.3f, beginning run for %d steps..\n", N, temperature, nsteps);

	for (t=0; t<nsteps; t++){
//		printf("Integrating\n");
		integrateLangevin(DT, temperature);

		if (calcNeigh()) { 
			rebuildCount ++;
//			printNeighCount(neighCount); 
#ifdef FLUSH			
			fflush(neighCount);
#endif
		}

		if (t% WFREQ == 0) {
//		if (1) {
			printf("step %d with %d rebuilds so far.\r",t,rebuildCount);

			printBubble(bubbles);

#ifdef FLUSH	
			fflush(stdout);	fflush(bubbles); 
#endif

			// print energy
			fprintf(energy, "%d\t%f\t%f\t%f\t%f\t%f\n",t, calcTemp(), intraE, interE, dihedralE, hardE );
			writeVTF(traj);
		}
		if (t% 1000000 == 0) writeRestart("restart");
	}
	
	printf("Done. Neighbour list was rebuilt about every %d steps.\n",nsteps/rebuildCount);

	free(seed);
	fclose(traj);
	fclose(minim);
	fclose(energy);
	fclose(neighCount);
	fclose(bubbles);

}

static float ziggurat(int thread_num) {  
	
	uint32_t  tmp = seed[thread_num];
      	float random = r4_nor (& tmp , kn, fn, wn );
	seed[thread_num] = tmp;
	return random;
      }

static void printMat(float matrix[][3]) {
	int i;
	for (i=0; i<2*N; i++)
		printf("%.2f, %.2f, %.2f\n", matrix[i][0], matrix[i][1], matrix[i][2]);
	printf("\n");
}

static void printBubble (FILE * bubbleFile) {
	int i;
	for (i=0; i<N; i++) 
		fprintf(bubbleFile, "%d ", isBound[i]) ;
	fprintf(bubbleFile, "\n");
}

static void zero(float matrix[][3]) {
	int i,j; 
	for (i=0;i<2*N;i++) {
		matrix[i][0]=0;
		matrix[i][1]=0;
		matrix[i][2]=0;
	}
}

static float calcTemp () {
	float kinE=0;
	int i;

	#pragma omp parallel for reduction(+:kinE)
	for (i=0; i<2*N; i++) {
		kinE += v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
	}

	return MASS*kinE/(3*2*N);
}

static float maxForce() {
	int i;
	float tmp, max = 0;

	for (i=0; i<2*N; i++) {
		tmp = f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];
		if (tmp>max) max=tmp;
	}
	return sqrt(max);
}
