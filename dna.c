#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "stdbool.h"
#include "omp.h"

#include "ziggurat_openmp.c"

#define N 1000
#define NSTEPS 10000000

#define INTRA_BOND_LENGTH 0.5
#define INTER_BOND_LENGTH 1.1

#define SKIN 0.1

#define INTER_BOND_CUT 1.2
#define HARD_CUT 0.4

#define K_BOND 100

#define PHI_1 -100
#define PHI_2 -95

#define K_DIHEDRAL 10
#define EPSILON_DIHEDRAL 1

#define MAX_NEIGH 2*N-1
#define CUT_NEIGH 2.0

#define GAMMA 0.1
#define TEMP 0.2

#define PI 3.1415

#define MASS K_BOND

static void printMat(float [][3]) ;
static void printNeigh();
static void zero(float [][3]) ;
static void genLadder();
static void genDNA(float);
static void genVel();
static float randn();
static float calcTemp();
static FILE * initVTF();
static void writeVTF(FILE *);
static void integrateNVE(float, int, int) ;
static void integrateLangevin(float, float ) ;
static void calcNeigh();
static float ziggurat(int num_thread) ;  


// global variables

float potE, kinE, temp;

float x[2*N][3];
float v[2*N][3];
float f[2*N][3];
int neigh[2*N][MAX_NEIGH+1];
bool isBound[N];

// used for ziggurat
static float fn[128];
static uint32_t kn[128];
static float wn[128];
// used for SHR3
static uint32_t *seed;

#include "dihedral.c"
#include "harmonic.c"
#include "hardcore.c"
#include "harcos.c"

static int rigid = 0;
#include "langevin.c"


void main() {
       
	// seeding a random stream for each thread
	r4_nor_setup ( kn, fn, wn );
	uint32_t  jsr = 123456;
  	
	int max_threads = omp_get_max_threads();
	seed = ( uint32_t * ) malloc ( max_threads * sizeof ( uint32_t ) );
//	printf("max threads : %d\n",max_threads);
//	omp_set_num_threads(1);
	
	for (int thread = 0; thread < max_threads; thread++ ) 
		seed[thread] = shr3_seeded ( &jsr );

//	genVel();
//	genLadder();
	genDNA(10.5);
	zero(f);
	zero(v);

	isBound[0]=1;
	isBound[N-1]=1;

	FILE *vtf = initVTF();
	FILE *bubbles = fopen("/tmp/bubbles.dat", "w");

	// minimization via Langevin at 0 temperature
	
	for (int t=0; t<100000; t++){
		integrateLangevin(0.001,0);
		if (t%1000 ==0)
			writeVTF(vtf);
	}

	for (int t=0; t<NSTEPS; t++){
//		printf("Integrating\n");
		integrateLangevin(0.1, TEMP);
//		integrateNVE(0.001,0,1);

		if (t%100 == 0) {
//			printf("Calculating neighbors\n");
			calcNeigh();
		}

		if (t%1000 ==0) {
//		if (1) {
	//		temp=calcTemp();
	//		printf("step %d, temp %f\n",t,temp);
			printf("\rstep %d",t);
			fflush(stdout);

			for (int i=0; i<N; i++) fprintf(bubbles,"%d ", isBound[i]) ;
			fprintf(bubbles, "\n");
			fflush(bubbles);

//			for (int i=0; i<N; i++) 
//				printf("%d ", neigh[i][0]);
//			printf("\n");

			
		//	printNeigh();
			writeVTF(vtf);
		}
	}
	
	printf("\n");

	free(seed);
	fclose(vtf);
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

static void zero(float matrix[][3]) {
	int i,j; 
	for (i=0;i<2*N;i++) {
		matrix[i][0]=0;
		matrix[i][1]=0;
		matrix[i][2]=0;
	}
}

/********************* Initial config functions *******************/

// generate a fully extended dna polymer

static void genLadder () {
	float yshift = -(N-1)*INTRA_BOND_LENGTH/2;

	for (int i=0; i<N; i++){

		x[2*i][0] = - INTER_BOND_LENGTH/2;
		x[2*i][1] = yshift + i* INTRA_BOND_LENGTH;
		x[2*i][2] = 0;
		
		x[2*i+1][0] = INTER_BOND_LENGTH/2;
		x[2*i+1][1] = yshift + i* INTRA_BOND_LENGTH;
		x[2*i+1][2] = 0;
	}
}

// generate a twisted DNA polymer

static void genDNA (float pitch) {
	float theta = 2*PI/pitch;
	float stepAngle = 0;

	// shift in y to center 
	float yshift = -(N-1)*(INTRA_BOND_LENGTH-0.1)/2;
	

	for (int i=0; i<N; i++){

		x[2*i][0] = INTER_BOND_LENGTH/2* cos(stepAngle+PI);
		x[2*i][1] = yshift + i* (INTRA_BOND_LENGTH-0.1);
		x[2*i][2] = INTER_BOND_LENGTH/2* sin(stepAngle+PI) ;
		
		x[2*i+1][0] = INTER_BOND_LENGTH/2 *cos(stepAngle) ;
		x[2*i+1][1] = yshift + i* (INTRA_BOND_LENGTH-0.1);
		x[2*i+1][2] = INTER_BOND_LENGTH/2 *sin (stepAngle);

		stepAngle+=theta;
	}
}
/********************** VTF functions *******************/

// initialize vtf file and return the pointer to the file

static FILE * initVTF() {
	FILE *vtf;

	vtf=fopen("/tmp/deneme.vtf", "w");

//	fprintf(vtf,"atom 0:%d radius %f name DNA\n", 2*N-1, HARD_CUT/2);

	//write atoms
	for (int i=0; i<2*N; i++)  
		fprintf(vtf,"a %d r %f c %d resid %d\n", i, HARD_CUT/2, i%2, i/2);
	
	// write bonds
	for (int i = 0 ; i<N; i++){
		//intra bonds
		if ( i < N-1) {
			fprintf(vtf,"bond %d:%d\n", 2*i, 2*i+2);
			fprintf(vtf,"bond %d:%d\n", 2*i+1, 2*i+3);
		}
		//inter bonds
		fprintf(vtf,"bond %d:%d\n", 2*i,2*i+1);
	}
	return vtf;
}

// write current timestep to the vtf file 

static void writeVTF(FILE *vtf) {
	int i;
	fprintf(vtf,"timestep\n");
	for (i=0; i<2*N; i++)
	  fprintf(vtf,"%.3f %.3f %.3f\n",x[i][0],x[i][1],x[i][2]);
}

/********************* Neighbour functions *******************/

// Neighbour list builder

static void calcNeigh() {
	int i,j;
	float del[3], rsq;
	float cutsq = CUT_NEIGH*CUT_NEIGH;

	// zero neighbour counts	
	for (i=0;i<2*N;i++) neigh[i][0]=0;

	#pragma omp parallel for private (j,del,rsq)
	for (i=0 ; i < 2*N; i++) for (j=i+1; j < 2*N; j++) {

		del[0]=x[i][0]-x[j][0];
		del[1]=x[i][1]-x[j][1];
		del[2]=x[i][2]-x[j][2];

		rsq= del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

		if (rsq<cutsq) 	neigh[i][ ++ neigh[i][0] ]=j;	
	}
}

// Neighbour list printer for debugging 

static void printNeigh() {
	int numNeigh;
	for (int i=0; i<2*N; i++) {
		numNeigh = neigh[i][0];
		printf("%d (%d) : ", i,numNeigh );
		for (int j = 1; j<=numNeigh; j++)
			printf("%d ",neigh[i][j]);
		printf("\n");
	}
	printf("\n");
}


/********************* Integrator functions *******************/

// Velocity Verlet integration

static void integrateNVE(float dt, int begin, int shift) {
	
	int j;

	float mult=0.5*dt/MASS;
	
	for (int i=begin; i<2*N; i+=shift) for (j=0; j<3; j++) {
		x[i][j]+= v[i][j]*dt + f[i][j]*mult*dt;
		v[i][j]+= f[i][j]*mult;
	}

//	calcForce();

	for (int i=begin; i<2*N; i+=shift) 
		for (j=0; j<3; j++) 
			v[i][j] += f[i][j]*mult;
}




static float calcTemp () {
	kinE=0;

	for (int i=0; i<2*N; i++) {
		kinE += v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
	}

	return MASS*kinE/(3*2*N);
}


