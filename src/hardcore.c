// hardcore repulsive potentials for non-bonded interactions
// based on mie n-m potentials

// mie(n,m,r) = c * k * ( (sigma/r)^n - (sigma/r)^m ) + k
// where 
// c     = (n/(n-m)) * (n/m)^(m/(n-m))
// sigma = r0*(m/n)**(1./(n-m))
// and
// if r>r0  	: hardcore_n_m(r) = mie(n,m,r)
// else 	: hardcore_n_m(r) = 0 

static const float hardCutSq = R_INTRA*R_INTRA;
static float sigma2_intra, sigma2_inter;
// hardcore repulsive potential based on mie 4-2 
// c = 4
// sigma(0.5) = 0.3535533905932738
// k = 0.05 * K_INTRA for best compatibility

static float hardcore_4_2(int i, int j, float sigma2) {

	float del[3], rsq;
	
	del[0]=x[i][0]-x[j][0];
	del[1]=x[i][1]-x[j][1];
	del[2]=x[i][2]-x[j][2];

	// this is squared distance
	rsq = del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

	// terminate early if beyond cutoff
	if ( rsq > hardCutSq ) return 0;

	float r2inv, fmult;

	r2inv  = 1/rsq;
	
	fmult  = 4.0 * r2inv * sigma2 * ( 4.0*sigma2*r2inv - 2.0) ;
	fmult *= K_INTRA * 0.05 * r2inv;

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

	return 4.0 * K_INTRA * 0.05 * sigma2 * r2inv * (sigma2 * r2inv - 1.0 ); 
}

// hardcore repulsive potential based on mie 4-2 
// c = 2.598
// sigma(0.5) = 0.37991784282579627
// k = 0.03 * K_INTRA for best compatibility

static float hardcore_6_2(int i, int j, float sigma2) {

	float del[3], rsq;
	
	del[0]=x[i][0]-x[j][0];
	del[1]=x[i][1]-x[j][1];
	del[2]=x[i][2]-x[j][2];

	// this is squared distance
	rsq = del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

	// terminate early if beyond cutoff
	if ( rsq > hardCutSq ) return 0;

	float r2inv, r4inv, fmult;
	float sigma4 = sigma2*sigma2;

	r2inv = 1/rsq;
	r4inv = r2inv*r2inv;
	
	fmult  = 2.598 * r2inv * sigma2 * ( 6.0*sigma4*r4inv - 2.0 ) ;
	fmult *= K_INTRA * 0.03 * r2inv;

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

	return 2.598 * K_INTRA * 0.03 * r2inv * sigma2 * (sigma4*r4inv - 1 ); 
}

// hardcore repulsive potential based on 8-4 mie potential
// c = 4.0
// sigma(0.5) = 0.42044820762685725
// k = 0.001 * K_INTRA

static float hardcore_8_4(int i, int j, float sigma2) {

	float del[3], rsq;
	
	del[0]=x[i][0]-x[j][0];
	del[1]=x[i][1]-x[j][1];
	del[2]=x[i][2]-x[j][2];

	// this is squared distance
	rsq = del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

	// terminate early if beyond cutoff
	if ( rsq > hardCutSq ) return 0;

	float r2inv, r4inv, fmult;
	float sigma4 = sigma2*sigma2;

	r2inv  = 1/rsq;
	r4inv  = r2inv*r2inv;
	
	fmult  = r4inv * sigma4 * 16 * (2*sigma4*r4inv - 1) ;
	fmult *=  0.001 * K_INTRA * r2inv;

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

	return 4 * 0.001 * K_INTRA * sigma4*r4inv*(sigma4*r4inv-1 ); 
}
// hardcore repulsive potential based on 12-6 mie potential
// c = 4.0
// sigma(0.5) = 0.44544935907016964
// k = 0.0036 * K_INTRA

static float hardcore_12_6(int i, int j, float sigma2) {

	float del[3], rsq;
	
	del[0]=x[i][0]-x[j][0];
	del[1]=x[i][1]-x[j][1];
	del[2]=x[i][2]-x[j][2];

	// this is squared distance
	rsq = del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

	// terminate early if beyond cutoff
	if ( rsq > hardCutSq ) return 0;

	float r2inv, r6inv, fmult;
	float sigma6 = sigma2*sigma2*sigma2;

	r2inv  = 1/rsq;
	r6inv  = r2inv*r2inv*r2inv;
	
	fmult  = r6inv * sigma6 * 24 * (2*sigma6*r6inv - 1) ;
	fmult *=  0.0036 * K_INTRA * r2inv;

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

	return 4 * 0.0036 * K_INTRA * sigma6*r6inv*(sigma6*r6inv-1 ); 
}


//quadratic softcore repulsive potential
// if r<rCut
// 	V(r) = K*(r-r0)^2
// 	F(r) = -dV/dr = -2*K*(r-r0)
// else
// 	V = 0 , F = 0

// squared cutoff used for performance

static float softcore(int i, int j, float r0) {

	float del[3], r;
	
	del[0]=x[i][0]-x[j][0];
	del[1]=x[i][1]-x[j][1];
	del[2]=x[i][2]-x[j][2];

	// this is squared distance
	r = del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

	// terminate early if beyond cutoff
	if ( r > r0 * r0 ) return 0;

	float dr,kdr,fmult;

	r  = sqrt(r);
	dr = r - r0;
	kdr= K_INTRA * dr;
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

