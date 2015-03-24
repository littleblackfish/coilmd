#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "signal.h"
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

// integrator parameters

#define DT 0.1
#define GAMMA 1.0
#define MASS 200

// intra strand parameters
#define R_INTRA 0.5
#define K_INTRA 200

// inter-strand parameters
#define R_INTER 1.1
#define K_INTER 100
#define E_INTER 1.0
#define CUT_INTER 1.2

// dihedral parameters
#define K_DIHED 1.0
#define E_DIHED 0.1
   
// non-bonded parameters
#define K_HARDCORE 400.0
#define CUT_NEIGH 1.0
#define MAX_NEIGH 64

// global energy variables
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

int N,N2,t,stopRunning=0;


// ladder case does not have any twist
#ifdef LADDER
static const float PHI_1 = 180.0;
static const float PHI_2 = 180.0;
#else
static const float PHI_1 = -100.0;
static const float PHI_2 = -95.0; 
#endif

static float sin1,cos1,sin2,cos2;

#include "helper.c"
#include "generators.c"
#include "restart.c"
#include "vtf.c"
#include "neighbour.c"
#include "langevin.c"

void graceful_exit(int signo)
{
	if (signo == SIGINT) printf("\nReceived SIGINT ");
	if (signo == SIGUSR1) printf("\nReceived SIGUSR1 ");
	if (signo == SIGUSR2) printf("\nReceived SIGUSR2 ");
	 
	printf ("at timestep %d, will abort after this step.",t);
	fflush(stdout);
	stopRunning = 1;
}

void main(int argc, char ** argv ) {
	signal(SIGINT,  graceful_exit);
	signal(SIGUSR1, graceful_exit);
	signal(SIGUSR2, graceful_exit);

	printf("MD simulation of a coarse grained DNA model.\nCompiled in ");
	#ifdef CIRCULAR 
	printf("CIRCULAR ");
	#endif
	#ifdef LADDER 
	printf("LADDER ");
	#else
	printf("COIL ");
	#endif
	printf("mode ");
	#ifdef _OPENMP
	printf("with OpenMP, running %d threads.\n", omp_get_max_threads());
	#else
	printf("without OpenMP.\n");
	#endif
	
	if (argc<4)  { printf("I need 3 parameters: N, nsteps, temperature.\n"); exit(1); }

	// get the parameters :  N, NSTEPS, TEMPERATURE
	N = atoi(argv[1]);
	N2 = 2*N;

	int nsteps = atoi(argv[2]);
	float temperature = atof (argv[3]);
	
	int i, rebuildCount = 0, rebuildDelay = 101; 

	// seeding a random stream for each thread
	r4_nor_setup ( kn, fn, wn );
	uint32_t  jsr = 123456;
  	
	int max_threads = omp_get_max_threads();
	seed = malloc ( max_threads * sizeof ( uint32_t ) );
	
	for (i=0; i < max_threads; i++ ) 
		seed[i] = shr3_seeded ( &jsr );


	// the order of hardcore repulsion and associated sigma
	float n=4,m=2;
	float sigma = R_INTRA * pow((m/n), (1.0/(n-m)));
	printf ("sigma is %f\n", sigma);
	sigma2 = pow(sigma, 2);
	sigma4 = pow(sigma, 4);
	sigma6 = pow(sigma, 6);

	// calculate sin and cos shifts for dihedrals once 
	sin1 = sin(PHI_1/180.*M_PI);
	cos1 = cos(PHI_1/180.*M_PI);
	sin2 = sin(PHI_2/180.*M_PI);
	cos2 = cos(PHI_2/180.*M_PI);

	//allocate global arrays 
	
	x = malloc(sizeof((*x))*2*N);
	v = malloc(sizeof((*v))*2*N);
	f = malloc(sizeof((*f))*2*N);
	
	xRef	= malloc(sizeof((*xRef))*2*N);
	neigh	= malloc(sizeof((*neigh))*2*N);
	isBound	= malloc(sizeof(int)*N);

	// force neighbour rebuild at the first step
	
	zero(xRef);
	xRef[0][0] = CUT_NEIGH;
	
#ifndef CIRCULAR	
	// first and last beads are initially (and always) bound 
	isBound[0]=1;
	isBound[N-1]=1;
#endif

	
	// if there is no restart file, 
	// do minimization via langevin at 0 temperature
	if ( !readRestart("restart") ) {
		FILE *minim = initVTF("minim.vtf"); 

		#if defined LADDER
		genLadder();
		#elif defined CIRCULAR
		genCircCoil(12.3);
		#else
		genCoil(12.3);
		#endif
		
		zero(f);
		zero(v);
		t=0;

		printf ("Minimizing...");

		do {
			if (t%1000 ==0) writeVTF(minim);	
			calcNeigh();
			integrateLangevin(0.005,0);
			t++;
		} while (maxForce(f) > 0.5 && t < 100000 );

		printf ("done after %d steps.\n",t);
		fclose(minim);
		
		//write minimized config to restart file
		writeRestart("restart");
	
	}

	// open files in append mode 

	FILE *traj   	= initVTF("traj.vtf");
	FILE *energy	= fopen("energy.dat", "a"); 
	FILE *bubbles	= fopen("bubbles.dat", "a");
	#ifdef NDEBUG			
	FILE *neighCount= fopen("neigh.dat", "a");
	#endif

	//force neighbour rebuild at first step
	xRef[0][0] = CUT_NEIGH;

	printf("N=%d, T=%.3f, beginning run for %d steps..\n", N, temperature, nsteps);
	fflush(stdout);

	/******* MAIN LOOP BEGINS HERE ********/

	for (t=0; t<nsteps; t++){
	
		// neighbour list rebuild with delay
		if ( calcNeigh() ) { 
			rebuildCount ++ ;
			#ifdef NDEBUG			
			printNeighCount(neighCount); fflush(neighCount);
			#endif
		}

		// integration
		integrateLangevin(DT, temperature);
			
		// stop running if necessary;
		if (stopRunning) break;

		// print bubble matrix
		if (t % 10000 == 0)   printBubble(bubbles);

		
		//write trajectory and energy
		if (t % 100000 == 0) {
			writeVTF(traj);
			printEnergy(energy);
		}
		
		//write restart file and flush buffers
		if (t % 1000000 == 0) {
			writeRestart("restart");
			fflush(bubbles); fflush(traj); fflush(energy);
			printf("\rWritten restart at step %d.",t); fflush(stdout);
		}
	}
	
	/******* MAIN LOOP ENDS HERE ******/

	//write final restart
	printBubble(bubbles);
	writeVTF(traj);
	printEnergy(energy);
	writeRestart("restart");
	
	printf("\nStopped at step %d.\nNeighbour list was rebuilt %d times.\n",t,rebuildCount);

	//close all files
	free(x); free(v); free(f); free(xRef); 
	free(isBound); free(neigh); free(seed);
	fclose(traj); fclose(energy); fclose(bubbles);
	#ifdef NDEBUG			
	fclose(neighCount);
	#endif

}
