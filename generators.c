/********************* Initial config functions *******************/

// generate a fully extended dna polymer

static void genLadder () {
	float yshift = -(N-1)*INTRA_BOND_LENGTH/2;
	int i;

	for (i=0; i<N; i++){

		x[2*i][0] = - INTER_BOND_LENGTH/2;
		x[2*i][1] = yshift + i* INTRA_BOND_LENGTH;
		x[2*i][2] = 0;
		
		x[2*i+1][0] = INTER_BOND_LENGTH/2;
		x[2*i+1][1] = yshift + i* INTRA_BOND_LENGTH;
		x[2*i+1][2] = 0;
	}
}

// generate a twisted DNA polymer

static void genDNA (float pitch) {
	float theta = 2*M_PI/pitch;
	float stepAngle = 0;

	// shift in y to center 
	float yshift = -(N-1)*(INTRA_BOND_LENGTH-0.1)/2;
	int i;
	

	for (i=0; i<N; i++){

		x[2*i][0] = INTER_BOND_LENGTH/2* cos(stepAngle+M_PI);
		x[2*i][1] = yshift + i* (INTRA_BOND_LENGTH-0.1);
		x[2*i][2] = INTER_BOND_LENGTH/2* sin(stepAngle+M_PI) ;
		
		x[2*i+1][0] = INTER_BOND_LENGTH/2 *cos(stepAngle) ;
		x[2*i+1][1] = yshift + i* (INTRA_BOND_LENGTH-0.1);
		x[2*i+1][2] = INTER_BOND_LENGTH/2 *sin (stepAngle);

		stepAngle+=theta;
	}
}

