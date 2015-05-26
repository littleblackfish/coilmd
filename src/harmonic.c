// harmonic potential 
// V(r) = K*(r-r0)^2
// F(r) = -dV/dr = -2*K*(r-r0)

static float harmonic(int i, int j, float k, float r0)
{
	float del[3], r, dr, kdr, fmult;
	
	del[0]=x[i][0]-x[j][0]; 
	del[1]=x[i][1]-x[j][1]; 
	del[2]=x[i][2]-x[j][2]; 

	//current bond length
	r = sqrt(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]);

	dr = r-r0;

	kdr = k*dr;
	
//	if (r <= 0.0) {
//		printf("r is 0 for harmonic between %d and %d. This is weird. \n",i,j);
//		return kdr*dr;
//	}

	fmult = -2.0*kdr/r;

	del[0]*=fmult;
	del[1]*=fmult;
	del[2]*=fmult;
	
	//apply forces
	
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

	//return energy
	
	return  kdr*dr;	
}
