static const float hardCutSq = R_INTRA*R_INTRA;

// hardcore repulsive potential based on 12-6 mie (lennard-jones) potential

static float sigma6;

static float hardcore(int i, int j) {

	float del[3], rsq;
	
	del[0]=x[i][0]-x[j][0];
	del[1]=x[i][1]-x[j][1];
	del[2]=x[i][2]-x[j][2];

	// this is squared distance
	rsq = del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

	// terminate early if beyond cutoff
	if ( rsq > hardCutSq ) return 0;

	float r2inv, r6inv, r12inv, fmult;

	r2inv  = 1/rsq;
	r6inv  = r2inv*r2inv*r2inv;
	
	fmult  = r6inv * sigma6 * 24 * (2*sigma6*r6inv - 1) ;
	fmult *= K_HARDCORE * r2inv;

	del[0]*=fmult;
	del[1]*=fmult;
	del[2]*=fmult;

	// apply forces
	
	#pragma omp atomic update
	f[i][0] += del[0];
	#pragma omp atomic update
	f[i][1] += del[1];
	#pragma omp atomic update
	f[i][2] += del[2];
	
	#pragma omp atomic update
	f[j][0] -= del[0];
	#pragma omp atomic update
	f[j][1] -= del[1];
	#pragma omp atomic update
	f[j][2] -= del[2];

	// return energy

	return 4 * K_HARDCORE *sigma6*r6inv*(sigma6*r6inv-1 ); 
}




//quadratic softcore repulsive potential
// if r<rCut
// 	V(r) = K*(r-r0)^2
// 	F(r) = -dV/dr = -2*K*(r-r0)
// else
// 	V = 0 , F = 0

// squared cutoff used for performance

static float softcore(int i, int j) {

	float del[3], r;
	
	del[0]=x[i][0]-x[j][0];
	del[1]=x[i][1]-x[j][1];
	del[2]=x[i][2]-x[j][2];

	// this is squared distance
	r = del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

	// terminate early if beyond cutoff
	if ( r > CUT_INTRA * CUT_INTRA ) return 0;

	float dr,kdr,fmult;

	r  = sqrt(r);
	dr = r - R_INTRA;
	kdr= K_HARDCORE * dr;
	fmult = -2.0*kdr/r;
	
	del[0]*=fmult;
	del[1]*=fmult;
	del[2]*=fmult;

	// apply forces
	
	#pragma omp atomic update
	f[i][0] += del[0];
	#pragma omp atomic update
	f[i][1] += del[1];
	#pragma omp atomic update
	f[i][2] += del[2];
	
	#pragma omp atomic update
	f[j][0] -= del[0];
	#pragma omp atomic update
	f[j][1] -= del[1];
	#pragma omp atomic update
	f[j][2] -= del[2];

	// return energy

	return kdr*dr; 
}

