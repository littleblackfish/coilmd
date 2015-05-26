/********************* Initial config functions *******************/

// norm of 3d vector

float norm (float x[3]) {return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);}


// rotate vector x around unit vector u by theta radians

static void rotate( float x[3], float u[3], float theta) {

	float c = cos(theta);
	float s = sin(theta);
	float rx = x[0], ry = x[1], rz = x[2];
	float ux = u[0], uy = u[1], uz = u[2];

	x[0] = ( ux*ux*(1-c) + c    ) *rx\
	     + ( ux*uy*(1-c) - uz*s ) *ry\
	     + ( ux*uz*(1-c) + uy*s ) *rz;

	x[1] = ( ux*uy*(1-c) + uz*s ) *rx\
	     + ( uy*uy*(1-c) + c    ) *ry\
	     + ( uy*uz*(1-c) - ux*s ) *rz;

	x[2] = ( uz*ux*(1-c) - uy*s ) *rx \
	     + ( uz*uy*(1-c) + ux*s ) *ry \
	     + ( uz*uz*(1-c) + c    ) *rz;
	
//	printf("%f\n", norm(x)); 
}

// generate a fully extended ladder

static void genLadder () {
	float yshift = -(N-1)*R_INTRA/2;
	int i;

	for (i=0; i<N; i++){

		x[2*i][0] = - R_INTER/2;
		x[2*i][1] = yshift + i* R_INTRA;
		x[2*i][2] = 0;
		
		x[2*i+1][0] = R_INTER/2;
		x[2*i+1][1] = yshift + i* R_INTRA;
		x[2*i+1][2] = 0;
	}
}

// generate a twisted ladder

static void genCoil (float pitch) {
	float theta = 2*M_PI/pitch;
	float stepAngle = 0;

	// shift in y to center 
	float yshift = -(N-1)*(R_INTRA-0.1)/2;
	int i;
	

	for (i=0; i<N; i++){

		x[2*i][0] = R_INTER/2* cos(stepAngle+M_PI);
		x[2*i][1] = yshift + i* (R_INTRA-0.1);
		x[2*i][2] = R_INTER/2* sin(stepAngle+M_PI) ;
		
		x[2*i+1][0] = R_INTER/2 *cos(stepAngle) ;
		x[2*i+1][1] = yshift + i* (R_INTRA-0.1);
		x[2*i+1][2] = R_INTER/2 *sin (stepAngle);

		stepAngle+=theta;
	}
}

// generate a circular twisted ladder

static void genCircCoil(float pitch) {



	float phi = 2*M_PI/N;
	float r   = 0.5*N*(R_INTRA-0.1)/M_PI;

	float theta = 2*M_PI/pitch;
	int i, j;
	float tangent[3];
	tangent[2]=0;

	zero(x);
	
	// place initial pair
	
	x[0][2] = R_INTER/2;
	x[1][2] =-R_INTER/2; 

	// twist following pairs 
	
	for (i=1; i<N; i++) {
		
		// copy last pair
		x[2*i][0] = x[2*i-2][0];
		x[2*i][1] = x[2*i-2][1];
		x[2*i][2] = x[2*i-2][2];

		x[2*i+1][0] = x[2*i-1][0];
		x[2*i+1][1] = x[2*i-1][1];
		x[2*i+1][2] = x[2*i-1][2];

		// tangent of the circle
		tangent[0] = cos(i*phi);
		tangent[1] = -sin(i*phi);

		// rotate aroung the tangent 
		rotate(x[2*i], tangent, -theta);
		rotate(x[2*i+1], tangent, -theta);
	}

	// displace all pairs to a circular path

	for (i=0; i<N; i++) {
		x[2*i][0]+= r*sin(i*phi);
		x[2*i][1]+= r*cos(i*phi);

		x[2*i+1][0]+= r*sin(i*phi);
		x[2*i+1][1]+= r*cos(i*phi);
	}

}
	
