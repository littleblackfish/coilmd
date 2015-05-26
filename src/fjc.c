#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "omp.h"

//K_b = 1

#define N 1000
#define mass 10
#define BOND_LENGTH 1
#define K_BOND 100
#define K_HARD 100
#define MAX_NEIGH N/2
#define CUT_NEIGH 1
#define CUT_HARD 1.0

#define GAMMA 0.1
#define TEMP 0.1


float potE;
float kinE;

float x[N][3];
float v[N][3];
float f[N][3];
int neigh[N][MAX_NEIGH+1];


void printMat(float matrix[N][3]) {
	int i;
	for (i=0; i<N; i++)
		printf("%.2f, %.2f, %.2f\n", matrix[i][0], matrix[i][1], matrix[i][2]);
	printf("\n");
}

void zero(float matrix[N][3]) {
	int i,j; 
	for (i=0;i<N;i++) for (j=0;j<3;j++) matrix[i][j]=0.;
}

void genVel(float v[N][3]) {
	int i, j;

	//generate vel
	//using a uniform distribution for now

	for (i=0; i<N; i++) 
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

// generate a fully extended polymer

void genPos (float x[][3]) {
	int i;
	for (i=0; i<N; i++){
		x[i][0]=0.;
		x[i][1]=0.;
		x[i][2]=i*BOND_LENGTH;
	}
}

// Temperature calculation
// Kb is taken as 1

static float calcTemp () {
	kinE=0;
	

//	#pragma omp parallel for reduction (+:kinE)
	for (int i=0; i<N; i++) {
		kinE += v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
	}

	return mass*kinE/(3*N);
}


static void writeVTF(FILE *vtf) {
	int i;
	fprintf(vtf,"timestep\n");
	for (i=0; i<N; i++)
	  fprintf(vtf,"%.3f %.3f %.3f\n",x[i][0],x[i][1],x[i][2]);
}


static void calcForce(int begin, int shift) {
	int i,j,k;
	float del[3];
	float norm;
	zero(f);
	potE=0;

	//there are N-1 bonds
	
//	#pragma omp parallel for private (del,norm) 
	for (int i=begin; i<N-1; i+=shift) {
		
		del[0]=x[i][0]-x[i+1][0]; 
		del[1]=x[i][1]-x[i+1][1]; 
		del[2]=x[i][2]-x[i+1][2]; 

		norm = sqrt(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]);

		del[0]/=norm;
		del[1]/=norm;
		del[2]/=norm;

		norm-=BOND_LENGTH;

	//	potE += 0.5*norm*norm*K_BOND;

		norm*=K_BOND;

		del[0]*=norm;
		del[1]*=norm;
		del[2]*=norm;

		#pragma omp atomic update
		f[i][0] -= del[0];
		#pragma omp atomic update
		f[i][1] -= del[1];
		#pragma omp atomic update
		f[i][2] -= del[2];

		#pragma omp atomic update
		f[i+1][0] += del[0];
		#pragma omp atomic update
		f[i+1][1] += del[1];
		#pragma omp atomic update
		f[i+1][2] += del[2];

	}

	// non-bonded
	
	float rsq, cutsq = CUT_HARD*CUT_HARD;

	#pragma omp parallel for private (k,j,rsq,del,norm) 
	for (i=begin; i<N; i+=shift) for (k=1; k<neigh[i][0]+1; k++) {

		j=neigh[i][k];

		del[0]=x[i][0]-x[j][0];
		del[1]=x[i][1]-x[j][1];
		del[2]=x[i][2]-x[j][2];

		rsq = del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

	//	printf("%.2f %.2f %.2f %.2f\n",del[0],del[1],del[2],rsq);
		
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
			
			#pragma omp atomic update 
			f[i][0] -= del[0];
			#pragma omp atomic update 
			f[i][1] -= del[1];
			#pragma omp atomic update 
			f[i][2] -= del[2];
			
			
			#pragma omp atomic update 
			f[j][0] += del[0];
			#pragma omp atomic update 
			f[j][1] += del[1];
			#pragma omp atomic update 
			f[j][2] += del[2];
		}

	}
}

// Velocity verlet integration

static void VVintegrate(float dt, int begin, int shift) {
	
	int j;

	float mult=0.5*dt/mass;
	
	for (int i=begin; i<N; i+=shift) for (j=0; j<3; j++) {
		x[i][j]+= v[i][j]*dt + f[i][j]*mult*dt;
		v[i][j]+= f[i][j]*mult;
	}

	calcForce(begin,shift);

	for (int i=begin; i<N; i+=shift) 
		for (j=0; j<3; j++) 
			v[i][j] += f[i][j]*mult;
}

// very inefficient box muller transform

float randn() {
	float x,y,s;
	x=drand48()*2-1;
	y=drand48()*2-1;
	s= x*x+y*y;
	if(s<1) 
		return x *sqrt(-2*log(s)/s);
	else 
		return randn();
}


// applying random forces for langevin dynamics

void randForce(float dt, int begin, int shift)  {
	const float mult = sqrt(2*TEMP*GAMMA*mass/dt);
	for (int i=begin; i<N;i+=shift) {
		f[i][0]+=randn()*mult;
		f[i][1]+=randn()*mult;
		f[i][2]+=randn()*mult;
	}
}


static void langevinIntegrate(float dt, int begin, int shift) {
	
	int j;

	const float halfdtgamma = 0.5*GAMMA*dt;
	const float halfdtgammanorm = 1/(1+halfdtgamma);
	const float halfdtMass = 0.5*dt/mass;
	
	for (int i=begin; i<N; i+=shift) for (j=0; j<3; j++) {

		// updating velocity by half a step
		// v = v(t+0.5dt)
		
		v[i][j] -= halfdtgamma*v[i][j];
		v[i][j] += halfdtMass*f[i][j];

		// updating position by a full step
		// x = x(t+dt)
		
		x[i][j]+= v[i][j]*dt ;

	}

	// calculate forces for the next timestep
	// f = f(t+dt)
	
	calcForce(begin,shift);
	randForce(dt,begin,shift);

	// final update on velocity 
	// v = v(t+dt)
	
	for (int i=begin; i<N; i+=shift) for (j=0; j<3; j++) { 
			v[i][j] += f[i][j]*halfdtMass;
			v[i][j] *= halfdtgammanorm;
	}
}


// build neighbour list

static void calcNeigh(int begin, int skip) {
	int i,j;
	float del[3], rsq;
	float cutsq = CUT_NEIGH*CUT_NEIGH;
	
	for (i=0;i<N;i++) neigh[i][0]=0;

	#pragma omp parallel for private (del,rsq,j)
	for (int i=begin ; i<N; i+= skip) for (j=i+1; j<N; j++) {

		del[0]=x[i][0]-x[j][0];
		del[1]=x[i][1]-x[j][1];
		del[2]=x[i][2]-x[j][2];

		rsq= del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

		if (rsq<cutsq) {
			neigh[i][0]++;
			neigh[i][neigh[i][0]]=j;
		}
	}
}

void main() {
	int i;

	genVel(v);
	genPos(x);
	zero(f);

	FILE *vtf;

	vtf=fopen("/tmp/deneme.vtf", "w");
	fprintf(vtf,"atom 0:%d radius %f name DNA\n", N-1, CUT_HARD/2);
	fprintf(vtf,"bond 0::%d\n", N-1);
	
	float temp;

	int nsteps = 10000;
	
	int t; 

	{	
	//int nthreads =  omp_get_num_threads();
	//int thisthread = omp_get_thread_num();

	int nthreads =  1;
	int thisthread =0; 

	printf ("starting simulation\nnumthreads : %d, thisthread : %d\n",nthreads,thisthread);

	for (int t=0; t<nsteps; t++){
		
		
		langevinIntegrate(0.001,thisthread,nthreads);
		
		if (t%10 == 0){ 
			calcNeigh(thisthread,nthreads);
		}
		
		if (t%10 ==0) {
	//	if (1) {
//			temp=calcTemp();
//			printf("step %d, temp %f\n",t,temp);
			
			writeVTF(vtf);
		}


	}
	}
	
	printf("\n");

}

