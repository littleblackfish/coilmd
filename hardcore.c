//quadratic hardcore repulsive potential
// if r<rCut
// 	V(r) = K*(r-r0)^2
// 	F(r) = -dV/dr = -2*K*(r-r0)
// else
// 	V = 0 , F = 0

static float hardcore(float f[][3], int i, int j, float k, float r0 ) {

	float del[3], r, dr, kdr, rsq, fmult;
	
	del[0]=x[i][0]-x[j][0];
	del[1]=x[i][1]-x[j][1];
	del[2]=x[i][2]-x[j][2];

	rsq = del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

	if ( rsq<hardCutSq ) {

		r  = sqrt(rsq);
		dr = r-r0;
		kdr= k*dr;
	
/*		if (r <= 0.0) {
			printf("r is %f for hardcore between %d and %d\n",r,i,j);
			printNeigh();
			exit(1);
			return kdr*dr;
		}
*/
		fmult = -2.0*kdr/r;
		
		del[0]*=fmult;
		del[1]*=fmult;
		del[2]*=fmult;

		//apply forces
		
//		#pragma omp atomic update
		f[i][0] += del[0];
//		#pragma omp atomic update
		f[i][1] += del[1];
//		#pragma omp atomic update
		f[i][2] += del[2];
		
//		#pragma omp atomic update
		f[j][0] -= del[0];
//		#pragma omp atomic update
		f[j][1] -= del[1];
//		#pragma omp atomic update
		f[j][2] -= del[2];

		//return energy

		return kdr*dr; 
	}

	else 
		return 0.0;
}

