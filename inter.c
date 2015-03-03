// this is a monolithic inter-strand potential that includes
// a harmonic cosine potential between the beads
// and 2 dihedrals between the strands
// which are scaled and cut-off in sync with the harmonic cosine potential
//
// previously, this was implemented as two seperate potentials in harcos.c and dihedral.c
// this resulted in 3 function calls and many redundant allocations and calculations
// this monolithic version does away with most of those. 

#define E_BOND 1
#define K_BOND 100
#define BOND_CUT 1.2
#define r0 1.1
#define rCut BOND_CUT

static float sin_shift1;
static float cos_shift1;

static float sin_shift2;
static float cos_shift2;

static float inter(int i2, int i3) {

	float energy;
  	float vb2x,vb2y,vb2z,rsq;


	// middle bond btw i2 and i3
	
	vb2x = x[i2][0] - x[i3][0];
	vb2y = x[i2][1] - x[i3][1];
	vb2z = x[i2][2] - x[i3][2];
	
	rsq = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
	
	// terminate early if the bond is broken
	
	if (rsq > BOND_CUT*BOND_CUT) return 0.0;
	
	// continue otherwise
  		
	float fmult,r,dr,rk;
	float kMult, dkMult;
	float f1[3], f2[3], f3[3], f4[3];

	r = sqrt(rsq);
	dr = r - r0;
		
	
	if (r<r0)	//this is the harmonic part
	{  
		rk = K_BOND * dr;
		fmult = -2.0*rk/r;
		energy = rk*dr - E_BOND;

		// dihedral is not scaled
		kMult = 1;	dkMult = 0;

	}
	
	else		// this is the cosine part
	{ 
		rk = M_PI / (rCut-r0);
		energy = E_BOND * 0.5 * ( 1 - cos(dr*rk)) - E_BOND;
		fmult =  -E_BOND * 0.5 * rk * sin(dr*rk)/r;	
		
		// dihedral is also scaled
		kMult  = 0.5 * ( 1 + cos(dr*rk)) ;
	    	dkMult = -0.5*rk * sin(dr*rk)  ;

	}
	
	//printf(" %f %f %f\n",vb2x*fmult, vb2y*fmult, vb2z*fmult );

	
	// forces on i2 and i3
	 	
	f2[0] = vb2x*fmult;
	f2[1] = vb2y*fmult;
	f2[2] = vb2z*fmult;
	  
	f3[0] = -vb2x*fmult;
	f3[1] = -vb2y*fmult;
	f3[2] = -vb2z*fmult;
	
	//////// FIRST DIHEDRAL BEGINS HERE ///////////	
	
	// first and last bonds
	float vb1x,vb1y,vb1z,vb3x,vb3y,vb3z;

	// other variables
	float ax,ay,az,bx,by,bz,rasq,rbsq,rinv,ra2inv,rb2inv,rabinv;
	float df,df1,fg,hg,fga,hgb,gaa,gbb;
	float dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx,dthy,dthz;
	float c,s,p,sx2,sy2,sz2;
	int i1,i4;

	// first dihedral
	i1 = (i2+N2-2)%N2;
	i4 = (i3+2)%N2;
	
	// 1st bond
	
	vb1x = x[i1][0] - x[i2][0];
	vb1y = x[i1][1] - x[i2][1];
	vb1z = x[i1][2] - x[i2][2];
	
	// 3rd bond
	
	vb3x = x[i4][0] - x[i3][0];
	vb3y = x[i4][1] - x[i3][1];
	vb3z = x[i4][2] - x[i3][2];
	
	// c,s calculation
	
	ax = vb1y*vb2z - vb1z*vb2y;
	ay = vb1z*vb2x - vb1x*vb2z;
	az = vb1x*vb2y - vb1y*vb2x;
	bx = vb3y*vb2z - vb3z*vb2y;
	by = vb3z*vb2x - vb3x*vb2z;
	bz = vb3x*vb2y - vb3y*vb2x;
	
	rasq = ax*ax + ay*ay + az*az;
	rbsq = bx*bx + by*by + bz*bz;
	
	rinv = ra2inv = rb2inv = 0.0;
	if (r > 0) rinv = 1.0/r;
	if (rasq > 0) ra2inv = 1.0/rasq;
	if (rbsq > 0) rb2inv = 1.0/rbsq;
	rabinv = sqrt(ra2inv*rb2inv);
	
	
	// cos(phi)
	c = (ax*bx + ay*by + az*bz)*rabinv;
	
	// sin(phi)
	s = r*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);
	
	// error check
	
	if (c > 1.0) c = 1.0;
	if (c < -1.0) c = -1.0;
	
	// V(phi) = f(rMid) * p(phi)
	// f(r)   = kMult (piecewise)
	// p(phi) = k * [ 1 - cos (phi - shift) ] - vShift 
	
	p = 1 - c*cos_shift1 - s*sin_shift1;
	df1 =   s*cos_shift1 - c*sin_shift1;
	
	p = K_DIHED*p - E_DIHED;
	
	energy += kMult * p;
	
	// force calculation begins here
	
	fg = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z;
	hg = vb3x*vb2x + vb3y*vb2y + vb3z*vb2z;
	fga = fg*ra2inv*rinv;
	hgb = hg*rb2inv*rinv;
	gaa = -ra2inv*r;
	gbb = rb2inv*r;
	
	dtfx = gaa*ax;
	dtfy = gaa*ay;
	dtfz = gaa*az;
	dtgx = fga*ax - hgb*bx;
	dtgy = fga*ay - hgb*by;
	dtgz = fga*az - hgb*bz;
	dthx = gbb*bx;
	dthy = gbb*by;
	dthz = gbb*bz;
	
	df = -kMult * K_DIHED * df1;
	
	sx2 = df*dtgx;
	sy2 = df*dtgy;
	sz2 = df*dtgz;
	
	f1[0] = df*dtfx;
	f1[1] = df*dtfy;
	f1[2] = df*dtfz;
	
	f2[0] += sx2 - f1[0] - vb2x*rinv*dkMult*p;
	f2[1] += sy2 - f1[1] - vb2y*rinv*dkMult*p;
	f2[2] += sz2 - f1[2] - vb2z*rinv*dkMult*p;
	
	f4[0] = df*dthx;
	f4[1] = df*dthy;
	f4[2] = df*dthz;
	
	f3[0] += -sx2 - f4[0] + vb2x*rinv*dkMult*p;
	f3[1] += -sy2 - f4[1] + vb2y*rinv*dkMult*p;
	f3[2] += -sz2 - f4[2] + vb2z*rinv*dkMult*p;
	
	// apply force to first and last atoms
	
	#pragma omp atomic update
      	f[i1][0] += f1[0];
	#pragma omp atomic update
      	f[i1][1] += f1[1];
	#pragma omp atomic update
      	f[i1][2] += f1[2];
      
	#pragma omp atomic update
    	f[i4][0] += f4[0];
	#pragma omp atomic update
    	f[i4][1] += f4[1];
	#pragma omp atomic update
	f[i4][2] += f4[2];
	
	// second dihedral
	i1 = (i2+2)%N2;
	i4 = (i3+N2-2)%N2;
	
	// 1st bond
	
	vb1x = x[i1][0] - x[i2][0];
	vb1y = x[i1][1] - x[i2][1];
	vb1z = x[i1][2] - x[i2][2];
	
	// 3rd bond
	
	vb3x = x[i4][0] - x[i3][0];
	vb3y = x[i4][1] - x[i3][1];
	vb3z = x[i4][2] - x[i3][2];
	
	// c,s calculation
	
	ax = vb1y*vb2z - vb1z*vb2y;
	ay = vb1z*vb2x - vb1x*vb2z;
	az = vb1x*vb2y - vb1y*vb2x;
	bx = vb3y*vb2z - vb3z*vb2y;
	by = vb3z*vb2x - vb3x*vb2z;
	bz = vb3x*vb2y - vb3y*vb2x;
	
	rasq = ax*ax + ay*ay + az*az;
	rbsq = bx*bx + by*by + bz*bz;
	
	rinv = ra2inv = rb2inv = 0.0;
	if (r > 0) rinv = 1.0/r;
	if (rasq > 0) ra2inv = 1.0/rasq;
	if (rbsq > 0) rb2inv = 1.0/rbsq;
	rabinv = sqrt(ra2inv*rb2inv);
	
	
	// cos(phi)
	c = (ax*bx + ay*by + az*bz)*rabinv;
	
	// sin(phi)
	s = r*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);
	
	// error check
	
	if (c > 1.0) c = 1.0;
	if (c < -1.0) c = -1.0;
	
	// V(phi) = f(rMid) * p(phi)
	// f(r)   = kMult (piecewise)
	// p(phi) = k * [ 1 - cos (phi - shift) ] - vShift 
	
	p = 1 - c*cos_shift2 - s*sin_shift2;
	df1 =   s*cos_shift2 - c*sin_shift2;
	
	p = K_DIHED*p - E_DIHED;
	
	energy += kMult * p;
	
	// force calculation begins here
	
	fg = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z;
	hg = vb3x*vb2x + vb3y*vb2y + vb3z*vb2z;
	fga = fg*ra2inv*rinv;
	hgb = hg*rb2inv*rinv;
	gaa = -ra2inv*r;
	gbb = rb2inv*r;
	
	dtfx = gaa*ax;
	dtfy = gaa*ay;
	dtfz = gaa*az;
	dtgx = fga*ax - hgb*bx;
	dtgy = fga*ay - hgb*by;
	dtgz = fga*az - hgb*bz;
	dthx = gbb*bx;
	dthy = gbb*by;
	dthz = gbb*bz;
	
	df = -kMult * K_DIHED * df1;
	
	sx2 = df*dtgx;
	sy2 = df*dtgy;
	sz2 = df*dtgz;
	
	f1[0] = df*dtfx;
	f1[1] = df*dtfy;
	f1[2] = df*dtfz;
	
	f2[0] += sx2 - f1[0] - vb2x*rinv*dkMult*p;
	f2[1] += sy2 - f1[1] - vb2y*rinv*dkMult*p;
	f2[2] += sz2 - f1[2] - vb2z*rinv*dkMult*p;
	
	f4[0] = df*dthx;
	f4[1] = df*dthy;
	f4[2] = df*dthz;
	
	f3[0] += -sx2 - f4[0] + vb2x*rinv*dkMult*p;
	f3[1] += -sy2 - f4[1] + vb2y*rinv*dkMult*p;
	f3[2] += -sz2 - f4[2] + vb2z*rinv*dkMult*p;
	
	// apply force to first and last atoms
	
	#pragma omp atomic update
      	f[i1][0] += f1[0];
	#pragma omp atomic update
      	f[i1][1] += f1[1];
	#pragma omp atomic update
      	f[i1][2] += f1[2];
      
	#pragma omp atomic update
    	f[i4][0] += f4[0];
	#pragma omp atomic update
    	f[i4][1] += f4[1];
	#pragma omp atomic update
	f[i4][2] += f4[2];

	// apply accumulated forces
	
	#pragma omp atomic update     	 
      	f[i2][0] += f2[0];
	#pragma omp atomic update
      	f[i2][1] += f2[1];
	#pragma omp atomic update
      	f[i2][2] += f2[2];
    
	#pragma omp atomic update
    	f[i3][0] += f3[0];
	#pragma omp atomic update
    	f[i3][1] += f3[1];
	#pragma omp atomic update
    	f[i3][2] += f3[2];

	return energy;
}
