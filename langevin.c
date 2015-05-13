#include "dihedral.c"
#include "harmonic.c"
#include "hardcore.c"
#include "harcos.c"
#include "angle.c"
#include "ziggurat_openmp.c"


// used for ziggurat
static float fn[128];
static uint32_t kn[128];
static float wn[128];
// used for SHR3
#ifdef _OPENMP
static uint32_t *seed;
#else
static uint32_t seed;
#endif


#ifdef _OPENMP
// normal random number generator using ziggurat method
static float ziggurat(int thread_num) {  
	
	uint32_t  tmp = seed[thread_num];
      	float random = r4_nor (& tmp , kn, fn, wn );
	seed[thread_num] = tmp;
	return random;
      }
#endif

float rk = M_PI / ( CUT_INTER - R_INTER );

// Velocity Verlet integration with Langevin EoM

static void integrateLangevin(float dt, float temperature) 
{
	const float halfdtgamma = 0.5*GAMMA*dt;
	const float halfdtgammanorm = 1/(1+halfdtgamma);
	const float halfdtMass = 0.5*dt/MASS;
	const float randFmult = sqrt(2*temperature*GAMMA*MASS/dt);

	float c,s; 
	float kMult, dkMult;
	
	// reset energy
	
	intraE=0;
	interE=0;
	angleE=0;
	dihedralE=0;
	intraHardE=0;
	interHardE=0;

	#pragma omp parallel 
	{
	#ifdef _OPENMP	
	int thread_num = omp_get_thread_num();
	int num_threads = omp_get_num_threads();
	#endif

	int i,j,k;
	float del[3],norm,rsq;
	float inter;

	#pragma omp for	schedule(static)
	for (i=0; i<2*N; i++) 
	{
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
		
		del[0] = v[i][0]*dt;
		del[1] = v[i][1]*dt;
		del[2] = v[i][2]*dt;

		x[i][0]+= del[0] ;
		x[i][1]+= del[1] ;
		x[i][2]+= del[2] ;

		// updating displacement since last neighbour rebuild
		
		xRef[i][0] += del[0] ;
		xRef[i][1] += del[1] ;
		xRef[i][2] += del[2] ;

	}

	
	// calculate forces for the next timestep
	// f = f(t+dt)
	
/********************* BEGIN FORCE CALCULATION ***********************/

	// Initialize with random forces instead of zeros 
	
	#pragma omp for	schedule(static)
	for (i=0; i<2*N; i++) {
		#ifdef _OPENMP
		f[i][0]=ziggurat(thread_num)*randFmult;
		f[i][1]=ziggurat(thread_num)*randFmult;
		f[i][2]=ziggurat(thread_num)*randFmult;
		#else
		f[i][0] = r4_nor (&seed , kn, fn, wn )*randFmult;
		f[i][1] = r4_nor (&seed , kn, fn, wn )*randFmult;
		f[i][2] = r4_nor (&seed , kn, fn, wn )*randFmult;
		#endif

	}

//	float randF [2*N]; if (t>890 && t<910) getMagnitude(f, randF);

	// private force array
/*	
	float f[2*N][3];
	int myflag[2*N];
	for (i=0; i<2*N; i++) myflag[i]=0;
	zero(f);
*/
	// Calculate forces from intra-strand bonds

	#pragma omp for	schedule(static) reduction(+:intraE)

	#ifdef CIRCULAR
	for (i=0; i<2*N; i++)	// circular case is periodical
		intraE += harmonic(i, (i+2)%(2*N), K_INTRA, R_INTRA);
	#else
	for (i=0; i<2*N-2; i++) // linear case has 2 less bonds
		intraE += harmonic(i, i+2, K_INTRA, R_INTRA);
	#endif

//	float intraF [2*N];	if (t>890 && t<910) getMagnitude(f, intraF);

	// Calculate forces from inter-strand interaction
	
	#pragma omp for reduction(+:interE,dihedralE) schedule(static)
	for (i=0; i<N; i++) {
		j=2*i;
		k=j+1;
		
		#ifndef CIRCULAR	
		if ( i == 0 ) { 
			harmonic(j, k, K_INTER, R_INTER);
			angleE += angle( j+2, j, k, kMult, dkMult );
			angleE += angle( k+2, k, j, kMult, dkMult );
			continue;
		}
		
		if ( i == N-1 ) {
			harmonic(j, k, K_INTER, R_INTER);
			angleE += angle( j-2, j, k, kMult, dkMult );
			angleE += angle( k-2, k, j, kMult, dkMult );
			continue;
		}
		#endif

		// circular case is periodical 
		inter = harcos(j, k, &c, &s);

		if  ( inter != 0 ) 
		{ 
			interE += inter;
			isBound[i] = 1;

			kMult  = 0.5 * ( 1 + c ) ;
		    	dkMult = -0.5 * rk * s  ;

			dihedralE += dihedral ((j+N2-2)%N2, j, k, (k+2)%N2, sin1, cos1, kMult, dkMult);
			dihedralE += dihedral ((j+2)%N2, j, k, (k+N2-2)%N2, sin2, cos2, kMult, dkMult);
				
			angleE += angle( (j+N2-2)%N2, j, k, kMult, dkMult );
			angleE += angle( (j+2)%N2   , j, k, kMult, dkMult );
			angleE += angle( (k+2)%N2   , k, j, kMult, dkMult );
			angleE += angle( (k+N2-2)%N2, k, j, kMult, dkMult );
		}

		else	isBound[i]=0;
	}

//	float interF[2*N]; if (t>890 && t<910) getMagnitude(f, interF);

	// Calculate forces form hard-core repulsion
	
	#pragma omp for	reduction(+:hardE) schedule(static)
	for (i=0; i<2*N; i++) 
	{
		for (k=1; k<neigh[i][0]+1; k++) 
		{
		j = neigh[i][k];

		// intra strand hardcore repulsion
		if ( (i+j)%2 == 0 ) 	
		//	intraHardE += hardcore_4_2 (i, j, sigma2_intra);
			intraHardE += softcore (i, j, R_HC_INTRA);

		// inter strand hardcore repulsion
		else			
		//	interHardE += hardcore_4_2 (i, j, sigma2_inter);
			interHardE += softcore (i, j, R_HC_INTER);

		}
	}
/*
	if (t>890 && t<910) {
		float nonbF[2*N];
		getMagnitude(f, nonbF);
		float forces[2*N][4];
		float max[4]={0,0,0,0};

		printf("step %d\nrand\tintra\tinter\tnonbond\n", t);
		for (i=0; i<2*N; i++) {
			forces[i][3]=abs(nonbF[i]-interF[i]);
			forces[i][2]=abs(interF[i]-intraF[i]);
			forces[i][1]=abs(intraF[i]-randF[i]);
			forces[i][0]=randF[i];
			for (j=0; j<4; j++) if (forces[i][j]>max[j]) max[j]=forces[i][j];
			
			printf("%d\t%f\t%f\t%f\t%f\n",i, forces[i][0], forces[i][1], forces[i][2], forces[i][3]);
		}
		printf("MAX\n%f\t%f\t%f\t%f\n", max[0], max[1], max[2], max[3]);

	}
*/


/****************** END OF FORCE CALCULATION **************************/
	
	
	// final update on velocity 
	// v = v(t+dt)
	#pragma omp for schedule(static)	
	for (i=0; i<2*N; i++)  { 

			v[i][0] += f[i][0]*halfdtMass ;
			v[i][1] += f[i][1]*halfdtMass ;
			v[i][2] += f[i][2]*halfdtMass ;
			
			v[i][0] *=halfdtgammanorm;
                        v[i][1] *=halfdtgammanorm;
	                v[i][2] *=halfdtgammanorm;
	}
	}
}
