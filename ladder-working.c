#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "omp.h"

#include "ziggurat_openmp.c"

#define N 10000
#define NSTEPS 1000
#define MASS 10
#define INTRA_BOND_LENGTH 0.5
#define INTER_BOND_LENGTH 1.1
#define K_BOND 100
#define K_HARD 100
#define MAX_NEIGH 2*N-1
#define CUT_NEIGH 1.2
#define CUT_HARD 0.4

#define GAMMA 0.1
#define TEMP 0.1


static void printMat(float [][3]) ;
static void printNeigh();
static void zero(float [][3]) ;
static void genLadder();
static void genVel();
static float randn();
static float calcTemp();
static void writeVTF(FILE *);
static void integrateNVE(float, int, int) ;
static void integrateLangevin(float, int, int) ;
static float ziggurat(int num_thread) ;  


// global variables

float potE;
float kinE;

float x[2*N][3];
float v[2*N][3];
float f[2*N][3];
int neigh[2*N][MAX_NEIGH+1];

// used for SHR3 generator
static float fn[128];
static uint32_t kn[128];
static float wn[128];


static void calcForce() {
	zero(f);
//#pragma omp parallel 
	{
	int i,j,k;
	float del[3];
	float norm;
	potE=0;

//	intra bonds
	
	#pragma omp for
	for (i=0; i<2*N-2; i++) {
		
		del[0]=x[i][0]-x[i+2][0]; 
		del[1]=x[i][1]-x[i+2][1]; 
		del[2]=x[i][2]-x[i+2][2]; 

		norm = sqrt(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]);

		del[0]/=norm;
		del[1]/=norm;
		del[2]/=norm;

		norm-=INTRA_BOND_LENGTH;

	//	potE += 0.5*norm*norm*K_BOND;

		norm*=K_BOND;

		del[0]*=norm;
		del[1]*=norm;
		del[2]*=norm;

		f[i][0] -= del[0];
		f[i][1] -= del[1];
		f[i][2] -= del[2];

		f[i+2][0] += del[0];
		f[i+2][1] += del[1];
		f[i+2][2] += del[2];
	}

	//inter	bonds
#pragma omp for
	for (i=0; i<N; i++) {
		
		del[0]=x[2*i][0]-x[2*i+1][0]; 
		del[1]=x[2*i][1]-x[2*i+1][1]; 
		del[2]=x[2*i][2]-x[2*i+1][2]; 

		norm = sqrt(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]);

		del[0]/=norm;
		del[1]/=norm;
		del[2]/=norm;

		norm-=INTER_BOND_LENGTH;
//		printf("%f\t%f\n",INTRA_BOND_LENGTH,norm);

	//	potE += 0.5*norm*norm*K_BOND;

		norm*=K_BOND;

		del[0]*=norm;
		del[1]*=norm;
		del[2]*=norm;

		f[2*i][0] -= del[0];
		f[2*i][1] -= del[1];
		f[2*i][2] -= del[2];

		f[2*i+1][0] += del[0];
		f[2*i+1][1] += del[1];
		f[2*i+1][2] += del[2];

	}


	// non-bonded
	
	float rsq, cutsq = CUT_HARD*CUT_HARD;

#pragma omp for
	for (i=0; i<2*N; i++) for (k=1; k<neigh[i][0]+1; k++) {
//		printf("%d ", neigh[i][0]);

		j=neigh[i][k];

		del[0]=x[i][0]-x[j][0];
		del[1]=x[i][1]-x[j][1];
		del[2]=x[i][2]-x[j][2];

		rsq = del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

		if ( rsq<cutsq ) {
			norm=sqrt(rsq);
			
//			printf ("%d %d interacting with r=%.2f %.2f\n",i,j,rsq,norm);

			del[0]/=norm;
			del[1]/=norm;
			del[2]/=norm;

			norm -= CUT_HARD;
			norm *= K_HARD;

			del[0]*=norm;
			del[1]*=norm;
			del[2]*=norm;
			
			f[i][0] -= del[0];
			f[i][1] -= del[1];
			f[i][2] -= del[2];
			
			f[j][0] += del[0];
			f[j][1] += del[1];
			f[j][2] += del[2];
		
		}

	}
	}
}

// build neighbour list

static void calcNeigh() {
	int i,j;
	float del[3], rsq;
	float cutsq = CUT_NEIGH*CUT_NEIGH;
	
	for (i=0;i<2*N;i++) neigh[i][0]=0;

	for (i=0 ; i < 2*N; i++) for (j=i+1; j < 2*N; j++) {

		del[0]=x[i][0]-x[j][0];
		del[1]=x[i][1]-x[j][1];
		del[2]=x[i][2]-x[j][2];

		rsq= del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

		if (rsq<cutsq) 	neigh[i][ ++ neigh[i][0] ]=j;	
	}
}

static uint32_t *seed;

void main() {
       
	// seeding a random stream for each thread
	r4_nor_setup ( kn, fn, wn );
	uint32_t  jsr = 123456;
  	
	int max_threads = omp_get_max_threads();
	seed = ( uint32_t * ) malloc ( max_threads * sizeof ( uint32_t ) );
	
	for (int thread = 0; thread < max_threads; thread++ ) 
		seed[thread] = shr3_seeded ( &jsr );




//	genVel();
	genLadder();
	zero(f);
	zero(v);

	FILE *vtf;

	vtf=fopen("/tmp/deneme.vtf", "w");
	fprintf(vtf,"atom 0:%d radius %f name DNA\n", 2*N-1, CUT_HARD/2);
	
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
	
	float temp;

	printf("\n");
	
	for (int t=0; t<NSTEPS; t++){
		integrateLangevin(0.001,0,1);
//		integrateNVE(0.001,0,1);
//		printMat(f);

		if (t%100 == 0) calcNeigh();
		//calcNeigh();

		if (t%1000 ==0) {
//		if (t>53550) {
//		if (1) {
	//		temp=calcTemp();
	//		printf("step %d, temp %f\n",t,temp);
			printf("step %d\r",t);
			
		//	printNeigh();
			writeVTF(vtf);
		}


	}
	
	printf("\n");

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

//generate initial velocities

void genVel() {
	int i, j;

	//generate vel
	//using a uniform distribution for now

	for (i=0; i<2*N; i++) 
		for (j=0; j<3; j++) {
			v[i][j]= drand48();
		}

	//remove translation
	
	for (j=0; j<3; j++) {
		float sum=0;
		for (i=0; i<N; i++)
			sum += v[i][j];
		
		sum/=N;
		for (i=0; i<N; i++)
			v[i][j]-=sum;
		
		sum=0;
	}

	//remove rotation
}

// very inefficient box muller transform

static float randn() {
	float x,y,s;
	x=drand48()*2-1;
	y=drand48()*2-1;
	s= x*x+y*y;
	if(s<1) 
		return x *sqrt(-2*log(s)/s);
	else 
		return randn();
}

static float calcTemp () {
	kinE=0;

	for (int i=0; i<2*N; i++) {
		kinE += v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
	}

	return MASS*kinE/(3*2*N);
}

static void writeVTF(FILE *vtf) {
	int i;
	fprintf(vtf,"timestep\n");
	for (i=0; i<2*N; i++)
	  fprintf(vtf,"%.3f %.3f %.3f\n",x[i][0],x[i][1],x[i][2]);
}


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

// Velocity Verlet integration

static void integrateNVE(float dt, int begin, int shift) {
	
	int j;

	float mult=0.5*dt/MASS;
	
	for (int i=begin; i<2*N; i+=shift) for (j=0; j<3; j++) {
		x[i][j]+= v[i][j]*dt + f[i][j]*mult*dt;
		v[i][j]+= f[i][j]*mult;
	}

	calcForce();

	for (int i=begin; i<2*N; i+=shift) 
		for (j=0; j<3; j++) 
			v[i][j] += f[i][j]*mult;
}


// Velocity Verlet integration with Langevin EoM

static void integrateLangevin(float dt, int begin, int shift) {
	
	const float halfdtgamma = 0.5*GAMMA*dt;
	const float halfdtgammanorm = 1/(1+halfdtgamma);
	const float halfdtMass = 0.5*dt/MASS;

//	#pragma omp parallel 
	{
	
	int thread_num = omp_get_thread_num();
	int num_threads = omp_get_num_threads();
	int i;

	
	const float randFmult = sqrt(2*TEMP*GAMMA*MASS/dt);
	
//	#pragma omp parallel for
	for (i=0; i<2*N; i++) {

		// updating velocity by half a step
		// v = v(t+0.5dt)
		
		v[i][0] -= halfdtgamma*v[i][0];
		v[i][1] -= halfdtgamma*v[i][1];
		v[i][2] -= halfdtgamma*v[i][2];

		v[i][0] += halfdtMass*f[i][0];
		v[i][1] += halfdtMass*f[i][1];
		v[i][2] += halfdtMass*f[i][2];

		// updating position by a full step
		// x = x(t+dt)
		
		x[i][0]+= v[i][0]*dt ;
		x[i][1]+= v[i][1]*dt ;
		x[i][2]+= v[i][2]*dt ;

	}

	// calculate forces for the next timestep
	// f = f(t+dt)
	
//	#pragma omp master
	calcForce();
	

	// Introduce random forces 
	

//	#pragma omp parallel for
	for (int i=0; i<2*N; i++) {
		f[i][0]+=ziggurat(omp_get_thread_num())*randFmult;
		f[i][1]+=ziggurat(omp_get_thread_num())*randFmult;
		f[i][2]+=ziggurat(omp_get_thread_num())*randFmult;
		
//		f[i][0]+=randn()*randFmult;
//		f[i][1]+=randn()*randFmult;
//		f[i][2]+=randn()*randFmult;
	}
	

	// final update on velocity 
	// v = v(t+dt)
	
	#pragma omp parallel for
	for (i=0; i<2*N; i++)  { 
			v[i][0] += f[i][0]*halfdtMass;
			v[i][1] += f[i][1]*halfdtMass;
			v[i][2] += f[i][2]*halfdtMass;

			v[i][0] *= halfdtgammanorm;
			v[i][1] *= halfdtgammanorm;
			v[i][2] *= halfdtgammanorm;
	}
	}
}


