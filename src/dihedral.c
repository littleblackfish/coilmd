static float dihedral(int i1, int i2, int i3, int i4, float sin_shift, float cos_shift, float kMult, float dkMult ) {
	
	float vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
	float energy,f1[3],f2[3],f3[3],f4[3];
	float ax,ay,az,bx,by,bz,rasq,rbsq,rg,rginv,ra2inv,rb2inv,rabinv;
	float df,df1,fg,hg,fga,hgb,gaa,gbb;
	float dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx,dthy,dthz;
	float c,s,p,sx2,sy2,sz2;
//	float kMult, dkMult;
	
	// 1st bond
	
	vb1x = x[i1][0] - x[i2][0];
	vb1y = x[i1][1] - x[i2][1];
	vb1z = x[i1][2] - x[i2][2];
	
	// 2nd bond
	
	vb2x = x[i3][0] - x[i2][0];
	vb2y = x[i3][1] - x[i2][1];
	vb2z = x[i3][2] - x[i2][2];
	
	vb2xm = -vb2x;
	vb2ym = -vb2y;
	vb2zm = -vb2z;
	
	//calculate length of middle bond (vb2)     
	rg = sqrt(vb2x*vb2x + vb2y*vb2y + vb2z*vb2z) ;
/*	kMult = 1;
	dkMult = 0;
	
	float dr,rk;
	//scale the coefficient if bond is further than equilibrium distance
	if (rg >  R_INTER  && rg <  CUT_INTER ) {
	    	rk = M_PI / ( CUT_INTER - R_INTER ); 	//scaling coefficient (use rasq temporarily for efficiency)
	    	dr = rg -  R_INTER ; 		//rdistance from r0   (use rbsq temporarily for efficiency)
	      	kMult  = 0.5 * ( 1 + cos(dr*rk)) ;
	    	dkMult = -0.5 * rk * sin(dr*rk)  ;
	    	//printf("Scaled dihedral due to bond btw %d and %d with r=%g, kScaled=%g\n" ,i2, i3, rMid,kScaled);
	      	//printf("Scaled dihedral r=%.2f, kMult=%.2f, dkMult=%.2f\n" ,rg,kMult,dkMult);
	}
	
	  //skip this dihedral if bond is broken
	else if (rg >=  CUT_INTER  )  {
	      //printf("Broken bond btw %d and %d with r=  %g\n" ,i2, i3, rg);
	      //printf("Broken bond with r=%.2f\n" , rg);
	      	printf ("Why are you calculating broken dihedrals, dummy?\n");
	      	return 0;   
	}
*/	
	// 3rd bond
	
	vb3x = x[i4][0] - x[i3][0];
	vb3y = x[i4][1] - x[i3][1];
	vb3z = x[i4][2] - x[i3][2];
	
	// c,s calculation
	
	ax = vb1y*vb2zm - vb1z*vb2ym;
	ay = vb1z*vb2xm - vb1x*vb2zm;
	az = vb1x*vb2ym - vb1y*vb2xm;

	bx = vb3y*vb2zm - vb3z*vb2ym;
	by = vb3z*vb2xm - vb3x*vb2zm;
	bz = vb3x*vb2ym - vb3y*vb2xm;
	
	rasq = ax*ax + ay*ay + az*az;
	rbsq = bx*bx + by*by + bz*bz;
	
	rginv = ra2inv = rb2inv = 0.0;
	if (rg   > 1e-8) rginv = 1.0/rg;
	if (rasq > 1e-8) ra2inv = 1.0/rasq;
	if (rbsq > 1e-8) rb2inv = 1.0/rbsq;
	rabinv = sqrt(ra2inv*rb2inv);
	
	
	// cos(phi)
	c = (ax*bx + ay*by + az*bz)*rabinv;
	
	// sin(phi)
	s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);
	
	// error check
	/*
	    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
	        printf(str,"Dihedral tolerance problem. " ) ;
	    }
	*/
	if (c > 1.0) c = 1.0;
	if (c < -1.0) c = -1.0;
	
	// V(phi) = f(rMid) * p(phi)
	// f(r)   = kMult (piecewise)
	// p(phi) = k * [ 1 - cos (phi - shift) ] - E_DIHED 
	
	p = 1 - c*cos_shift - s*sin_shift;
	df1 =   s*cos_shift - c*sin_shift;
	
	p = K_DIHED * p - E_DIHED;
	
	energy = kMult * p;
	
	// force calculation begins here
	
	fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
	hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
	fga = fg*ra2inv*rginv;
	hgb = hg*rb2inv*rginv;
	gaa = -ra2inv*rg;
	gbb = rb2inv*rg;
	
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
	
	f2[0] = sx2 - f1[0] - vb2xm*rginv*dkMult*p;
	f2[1] = sy2 - f1[1] - vb2ym*rginv*dkMult*p;
	f2[2] = sz2 - f1[2] - vb2zm*rginv*dkMult*p;
	
	f4[0] = df*dthx;
	f4[1] = df*dthy;
	f4[2] = df*dthz;
	
	f3[0] = -sx2 - f4[0] + vb2xm*rginv*dkMult*p;
	f3[1] = -sy2 - f4[1] + vb2ym*rginv*dkMult*p;
	f3[2] = -sz2 - f4[2] + vb2zm*rginv*dkMult*p;
	
	// apply force to each of 4 atoms
	
	#pragma omp atomic update
      	f[i1][0] += f1[0];
	#pragma omp atomic update
      	f[i1][1] += f1[1];
	#pragma omp atomic update
      	f[i1][2] += f1[2];
      
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
    
	#pragma omp atomic update
    	f[i4][0] += f4[0];
	#pragma omp atomic update
    	f[i4][1] += f4[1];
	#pragma omp atomic update
	f[i4][2] += f4[2];
	
	return energy;
}
