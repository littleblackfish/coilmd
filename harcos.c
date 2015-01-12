
static float harcos(int i, int j, float k, float r0, float epsilon, float rCut) {	
	float energy;
  	float delx,dely,delz,fmult;
  	float rsq,r,dr,rk;
  	float mult;

	delx = x[i][0] - x[j][0];
	dely = x[i][1] - x[j][1];
	delz = x[i][2] - x[j][2];
	
	rsq = delx*delx + dely*dely + delz*delz;
	
	if (rsq < rCut*rCut) {

		r = sqrt(rsq);
	
		if (r<r0) {  		//this is the harmonic part
			dr = r - r0;
			rk = k * dr;
	
			fmult = -2.0*rk/r;
	
			energy = rk*dr - epsilon;
		}
	
		else if (r<rCut) { 	//this is the cosine part
			dr = r - r0;
	
			rk = PI / (rCut-r0);
	    
			energy = epsilon * 0.5 * ( 1 - cos(dr*rk)) - epsilon;
	
			fmult =  -epsilon * 0.5 * rk * sin(dr*rk)/r;
	
		}
	}
	
	else     {
//		printf("The interstrand bond %d is broken\n", i/2);
		return 0.0 ;
	}
	
	//printf(" %f %f %f\n",delx*fmult, dely*fmult, delz*fmult );
	
	// apply forces
	 	
	delx *= fmult;
	dely *= fmult;
	delz *= fmult;
	
	#pragma omp atomic update
	f[i][0] += delx;
	#pragma omp atomic update
	f[i][1] += dely;
	#pragma omp atomic update
	f[i][2] += delz;
	  
	#pragma omp atomic update
	f[j][0] -= delx;
	#pragma omp atomic update
	f[j][1] -= dely;
	#pragma omp atomic update
	f[j][2] -= delz;
	
	return energy;
	
}
