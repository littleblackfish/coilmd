static float angle(int i1, int i2 , int i3, float kMult, float dkMult)
{
	float delx1,dely1,delz1,delx2,dely2,delz2;
	float eangle,f1[3],f3[3];
	float dcostheta,tk;
	float rsq1,rsq2,r1,r2,c,a,a11,a12,a22;
	
	// 1st bond
	delx1 = x[i1][0] - x[i2][0];
	dely1 = x[i1][1] - x[i2][1];
	delz1 = x[i1][2] - x[i2][2];
	
	rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
	r1 = sqrt(rsq1);
	
	// 2nd bond
	
	delx2 = x[i3][0] - x[i2][0];
	dely2 = x[i3][1] - x[i2][1];
	delz2 = x[i3][2] - x[i2][2];
	
	rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
	r2 = sqrt(rsq2);
	
	// angle (cos and sin)
	
	c = delx1*delx2 + dely1*dely2 + delz1*delz2;
	c /= r1*r2;
	
	if (c > 1.0) c = 1.0;
	if (c < -1.0) c = -1.0;
	
	// force & energy
	
	dcostheta = c - cos(THETA0);
	tk = K_ANGLE * dcostheta;

	// unscaled energy
	
	eangle = tk*dcostheta - E_ANGLE;
	
	a = 2.0 * tk * kMult;

	a11 = a*c / rsq1;
	a12 = -a / (r1*r2);
	a22 = a*c / rsq2;
	
	f1[0] = a11*delx1 + a12*delx2;
	f1[1] = a11*dely1 + a12*dely2;
	f1[2] = a11*delz1 + a12*delz2;
	
	// last term is the correction for scaling
	
	f3[0] = a22*delx2 + a12*delx1 - dkMult*eangle*delx2/r2;
	f3[1] = a22*dely2 + a12*dely1 - dkMult*eangle*dely2/r2;
	f3[2] = a22*delz2 + a12*delz1 - dkMult*eangle*delz2/r2;
	
	// apply force to each of 3 atoms
	
	f[i1][0] += f1[0];
	f[i1][1] += f1[1];
	f[i1][2] += f1[2];
	
	f[i2][0] -= f1[0] + f3[0];
	f[i2][1] -= f1[1] + f3[1];
	f[i2][2] -= f1[2] + f3[2];
	
	f[i3][0] += f3[0];
	f[i3][1] += f3[1];
	f[i3][2] += f3[2];

	// return scaled energy 
	
	return eangle*kMult;
}
