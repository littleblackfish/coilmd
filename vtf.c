
// initialize vtf file and return the pointer to the file

static FILE * initVTF( char filename[]) {
	FILE *vtf;
	int i;

	// attempt to read vtf file
	vtf=fopen(filename, "r");

	// if file exists, append to it
	if (vtf) {
		fclose(vtf);
		vtf=fopen(filename, "a");
	}
	// if the file does not exist, initialize it
	else { 
		vtf = fopen(filename, "w");
	
		//write atoms
		for (i=0; i<2*N; i++)  
			fprintf(vtf,"a %d r %f c %d resid %d resname dna%d\n",i,R_INTRA/2, i%2, i/2, (i%2)+1);
		
		// write bonds
		for (i = 0 ; i<N; i++){
			//intra bonds
			if ( i < N-1) {
				fprintf(vtf,"bond %d:%d\n", 2*i, 2*i+2);
				fprintf(vtf,"bond %d:%d\n", 2*i+1, 2*i+3);
			}
			//inter bonds
			fprintf(vtf,"bond %d:%d\n", 2*i,2*i+1);
		}
	}
		
	return vtf;
}

// write current timestep to the vtf file 

static void writeVTF(FILE *vtf) {
	int i;
	fprintf(vtf,"timestep\n");
	for (i=0; i<2*N; i++)
	  fprintf(vtf,"%.3f %.3f %.3f\n",x[i][0],x[i][1],x[i][2]);
}


