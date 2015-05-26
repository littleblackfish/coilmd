static void printMat(float matrix[][3]) {
	int i;
	for (i=0; i<2*N; i++)
		printf("%.2f, %.2f, %.2f\n", matrix[i][0], matrix[i][1], matrix[i][2]);
	printf("\n");
}

static void getMagnitude(float f[][3], float mag[]) {
	int i;
	for (i=0; i<2*N; i++)
		mag[i] = sqrt(f[i][0]*f[i][0]+f[i][1]*f[i][1]+f[i][2]*f[i][2]);
}


static void printBubble (FILE * bubbleFile) {
	int i;
	for (i=0; i<N; i++) 
		fprintf(bubbleFile, "%d ", isBound[i]) ;
	fprintf(bubbleFile, "\n");
}


static void zero(float matrix[][3]) {
	int i,j; 
	
	for (i=0;i<2*N;i++) {
		matrix[i][0]=0;
		matrix[i][1]=0;
		matrix[i][2]=0;
	}
}

static float calcTemp (float v[][3]) {
	float kinE=0;

	int i;
	#pragma omp parallel for reduction(+:kinE) private(i)
	for ( i=0; i<2*N; i++) {
		kinE += v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
	}

	return MASS*kinE/(3*2*N);
}

static float maxForce(float f[][3]) {
	int i;
	float tmp, maxf = 0;

	#pragma omp parallel for reduction(max:maxf) private(tmp)
	for (i=0; i<2*N; i++) {
		tmp = f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];
		if (tmp>maxf) maxf=tmp;
	}
	return sqrt(maxf);
}

static void printEnergy (FILE * energy) {
			fprintf(energy, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n",t, calcTemp(v), intraE, interE, angleE, dihedralE, intraHardE, interHardE );
}
