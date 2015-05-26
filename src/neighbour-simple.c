// Neighbour list builder

static int calcNeigh() {
	int i,j;
	float del[3], rsq;

	// zero neighbour counts	
	for (i=0;i<2*N;i++) neigh[i][0]=0;

	#pragma omp parallel for private (j,del,rsq)
	for (i=0 ; i < 2*N; i++) for (j=i+1; j < 2*N; j++) {

		del[0]=x[i][0]-x[j][0];
		del[1]=x[i][1]-x[j][1];
		del[2]=x[i][2]-x[j][2];

		rsq= del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

		if (rsq<cutsq) 	neigh[i][ ++ neigh[i][0] ]=j;	
	}
}

// Neighbour list printer for debugging 

static void printNeigh() {
	int numNeigh;
	for (int i=0; i<2*N; i++) {
		numNeigh = neigh[i][0];
		printf("%d (%d) : ", i,numNeigh );
		for (int j = 1; j<=numNeigh; j++)
			printf("%d ",neigh[i][j]);
		printf("\n");
	}
	printf("\n");
}


