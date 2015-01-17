
// initialize vtf file and return the pointer to the file

static FILE * initVTF( char filename[]) {
	FILE *vtf;
	int i;

	vtf=fopen(filename, "w");

//	fprintf(vtf,"atom 0:%d radius %f name DNA\n", 2*N-1, HARD_CUT/2);

	//write atoms
	for (i=0; i<2*N; i++)  
		fprintf(vtf,"a %d r %f c %d resid %d\n", i, HARD_CUT/2, i%2, i/2);
	
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
	return vtf;
}

// write current timestep to the vtf file 

static void writeVTF(FILE *vtf) {
	int i;
	fprintf(vtf,"timestep\n");
	for (i=0; i<2*N; i++)
	  fprintf(vtf,"%.3f %.3f %.3f\n",x[i][0],x[i][1],x[i][2]);
}


