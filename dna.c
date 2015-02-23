#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#ifdef _OPENMP
  #include "omp.h"
#else 
   inline int omp_get_thread_num() {return 0;}
   inline int omp_get_num_threads() {return 1;}
   inline int omp_get_max_threads() {return 1;}
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// some parameters describing the model

#define DT 0.1
#define GAMMA 1

#define INTRA_BOND_LENGTH 0.5
#define INTER_BOND_LENGTH 1.1

#define INTER_BOND_CUT 1.2
#define HARD_CUT 0.4

#define K_BOND 100
#define MASS K_BOND

#define K_DIHED 1
#define E_DIHED 0.1

#define NEIGH_CUT 1.2
#define MAX_NEIGH 64

// energy variables
float intraE, interE, dihedralE, hardE;

// position, velocity, force 
float (*x)[3];
float (*v)[3];
float (*f)[3];

// displacement tracking for neighbour search
float (*xRef)[3];

// neighbour list
int (*neigh)[MAX_NEIGH+1];

// association data
int *isBound;

int N;

// ladder case does not have any twist
#ifdef LADDER
static const float PHI_1 = 180.0;
static const float PHI_2 = 180.0;
#else
static const float PHI_1 = -100.0;
static const float PHI_2 = -95.0; 
#endif

#include "helper.c"
#include "generators.c"
#include "restart.c"
#include "vtf.c"
#include "neighbour.c"
#include "langevin.c"

void main(int argc, char ** argv ) {

	#ifdef _OPENMP
		printf("Compiled with OpenMP, running %d threads.\n", omp_get_max_threads());
	#else
		printf("Compiled without OpenMP.\n");
	#endif
	
	if (argc<4)  {
		printf("Not enough parameters. Cannot run without enough parameters.\n");
		exit(1);
	}

	// get the parameters :  N, NSTEPS, TEMPERATURE
	N = atoi(argv[1]);
	int nsteps = atoi(argv[2]);
	float temperature = atof (argv[3]);
	
	int i, t=0, rebuildCount = 0, rebuildDelay = 101; 

	// seeding a random stream for each thread
	r4_nor_setup ( kn, fn, wn );
	uint32_t  jsr = 123456;
  	
	int max_threads = omp_get_max_threads();
	seed = malloc ( max_threads * sizeof ( uint32_t ) );
	
	for (i=0; i < max_threads; i++ ) 
		seed[i] = shr3_seeded ( &jsr );

	x = malloc(sizeof((*x))*2*N);
	v = malloc(sizeof((*v))*2*N);
	f = malloc(sizeof((*f))*2*N);
	
	xRef = malloc(sizeof((*xRef))*2*N);

	neigh = malloc(sizeof((*neigh))*2*N);

	isBound = malloc(sizeof(int)*N);

	isBound[0]=1;
	isBound[N-1]=1;

	FILE *traj  = initVTF("traj.vtf");
	FILE *energy  =	fopen("energy.dat", "a"); 
	FILE *bubbles =	fopen("bubbles.dat", "a");
#ifdef NDEBUG			
	FILE *neighCount   = fopen("neigh.dat", "a");
#endif
	
	// minimization via Langevin at 0 temperature
	if ( !readRestart("restart") ) {
		FILE *minim = initVTF("minim.vtf"); 

#if defined LADDER
		genLadder();
#elif defined CIRCULAR
		genCircCoil(12);
#else
		genCoil(12);
#endif
		
		zero(f);
		zero(v);

		printf ("Minimizing...");

		do {
			if (t%1000 ==0) writeVTF(minim);	
			calcNeigh();
			integrateLangevin(0.005,0);
			t++;
		} while (maxForce(f) > 0.5 && t < 100000 );

		printf ("done after %d steps.\n",t);
		fclose(minim);
		writeRestart("restart");
	}

	xRef[0][0]=1000; //force neighbour rebuild at first step

	printf("N=%d, T=%.3f, beginning run for %d steps..\n", N, temperature, nsteps);
	fflush(stdout);

	/******* MAIN LOOP BEGINS HERE ********/

	for (t=0; t<nsteps; t++){
	
		// neighbour list rebuild with delay
		if (rebuildDelay++ > 50 && calcNeigh()) { 
			rebuildCount ++ ;
			rebuildDelay = 0 ; 
#ifdef NDEBUG			
			printNeighCount(neighCount); 
			fflush(neighCount);
#endif
		}

		// integration
		integrateLangevin(DT, temperature);

		// print bubble matrix
		if (t % 10000 == 0)   printBubble(bubbles);

		
		//write trajectory and energy
		if (t % 100000 == 0) {
			writeVTF(traj);
			printEnergy(energy, t);
		}
		
		//write restart file and flush buffers
		if (t % 1000000 == 0) {
			writeRestart("restart");
			fflush(bubbles); fflush(traj); fflush(energy);
			printf("written restart at step %d.\r",t); fflush(stdout);
		}
	}
	
	/******* MAIN LOOP ENDS HERE ******/

	//write final restart
	writeRestart("restart");

	
	printf("Done. Neighbour list was rebuilt about every %d steps.\n",nsteps/rebuildCount);

	//close all files
	free(x);
	free(v);
	free(f);
	free(xRef);
	free(isBound);
	free(neigh);
	free(seed);
	fclose(traj);
	fclose(energy);
	fclose(bubbles);
#ifdef NDEBUG			
	fclose(neighCount);
#endif

}
