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


static void printMat(float [][3]) ;
static void printBubble (FILE *) ;
static void zero(float [][3]) ;
static float calcTemp();
static float maxForce();

// global variables

float intraE, interE, dihedralE, hardE;

float x[2*N][3];
float xRef[2*N][3];
float v[2*N][3];
float f[2*N][3];
int neigh[2*N][MAX_NEIGH+1];
int isBound[N];

#ifdef LADDER
static const float PHI_1 = 180.0;
static const float PHI_2 = 180.0;
#else
static const float PHI_1 = -100.0;
static const float PHI_2 = -95.0; 
#endif

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
	
	int i, t=0, rebuildCount = 0, rebuildDelay = 101; 

	// seeding a random stream for each thread
	r4_nor_setup ( kn, fn, wn );
	uint32_t  jsr = 123456;
  	
	int max_threads = omp_get_max_threads();
	seed = ( uint32_t * ) malloc ( max_threads * sizeof ( uint32_t ) );
	
	for (i=0; i < max_threads; i++ ) 
		seed[i] = shr3_seeded ( &jsr );


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
		} while (maxForce() > 0.1 && t <1000000 );

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

		// integration..
		integrateLangevin(DT, temperature);

		// print bubble matrix
		if (t % 10000 == 0)   printBubble(bubbles);

		//write restart file
		if (t % 1000000 == 0) {
			writeRestart("restart");
			fflush(bubbles);
//			fflush(traj);
//			fflush(energy);
		}
		
		//write trajectory and energy
		if (t % 100000 == 0) {
			writeVTF(traj);
			fprintf(energy, "%d\t%f\t%f\t%f\t%f\t%f\n",t, calcTemp(), intraE, interE, dihedralE, hardE );

#ifdef FLUSH	
			printf("step %d with %d rebuilds so far.\r",t,rebuildCount);
			fflush(stdout);	fflush(bubbles); fflush(energy);
#endif
		}
	}
	
	/******* MAIN LOOP ENDS HERE ******/

	//write final restart
	writeRestart("restart");

	
	printf("Done. Neighbour list was rebuilt about every %d steps.\n",nsteps/rebuildCount);

	//close all files
	free(seed);
	fclose(traj);
	fclose(energy);
	fclose(bubbles);
#ifdef NDEBUG			
	fclose(neighCount);
#endif

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
	#pragma omp parallel for reduction(+:kinE) private(i)
	for ( i=0; i<2*N; i++) {
		kinE += v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
	}

	return MASS*kinE/(3*2*N);
}

static float maxForce() {
	int i;
	float tmp, maxf = 0;

	#pragma omp parallel for reduction(max:maxf) private(tmp)
	for (i=0; i<2*N; i++) {
		tmp = f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];
		if (tmp>maxf) maxf=tmp;
	}
	return sqrt(maxf);
}
