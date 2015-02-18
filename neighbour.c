// Neighbour list builder
// xRef[N][3] is used to keep the positions from last rebuild. 
// everything closer than NEIGH_CUT are included

// some constants for performance 
static const float neighCutSq  = NEIGH_CUT * NEIGH_CUT ;
static const float hardCutSq = HARD_CUT*HARD_CUT;
static const float neighSkinSq = (NEIGH_CUT-HARD_CUT) * (NEIGH_CUT - HARD_CUT) ;

static int calcNeigh() {
	
	int rebuild = 0;

	#pragma omp parallel 
	{
	int i,j;
	float del[3], rsq;
	
	#pragma omp for reduction(||:rebuild) schedule(static)	
	for (i=0; i < 2*N; i++) {

		if (!rebuild) {

			rsq = xRef[i][0]*xRef[i][0] + xRef[i][1]*xRef[i][1] + xRef[i][2]*xRef[i][2];
			
			if (rsq > neighSkinSq) {
				rebuild+=1;
			}
		}
	}
	
	if (rebuild) {
		#pragma omp for schedule(guided,100) 
		for (i=0 ; i < 2*N; i++) {
			
			// reset reference 
			xRef[i][0] = 0.0;
			xRef[i][1] = 0.0;
			xRef[i][2] = 0.0;
			
			//zero neighbour counts
			neigh[i][0] = 0;

			//calculate new neighbors
			for (j=i+1; j < 2*N; j++) {

				del[0]=x[i][0]-x[j][0];
				del[1]=x[i][1]-x[j][1];
				del[2]=x[i][2]-x[j][2];

				rsq= del[0]*del[0]+del[1]*del[1]+del[2]*del[2];

				if (rsq<neighCutSq) {
					neigh[i][ ++ neigh[i][0] ] = j;	
				}
			}
		}
	}
	}

	return rebuild;
}

// Neighbour list printer for debugging 

static void printNeigh() {
	int numNeigh,i,j;
	for (i=0; i<2*N; i++) {
		numNeigh = neigh[i][0];
		printf("%d (%d) : ", i,numNeigh );
		for (j = 1; j<=numNeigh; j++)
			printf("%d ",neigh[i][j]);
		printf("\n");
	}
	printf("\n");
}

// Neighbour count printer for load balance benchmarking

static void printNeighCount(FILE *f) {
	int i;
	for (i=0; i<N; i++) 
		fprintf(f,"%d ", neigh[i][0]);
	fprintf(f,"\n");
}




