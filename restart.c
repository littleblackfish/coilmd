static void readRestart(char * filename) {
	int i,n;
	char check;
	FILE *restart = fopen(filename, "r");
	
	fscanf(restart,"%d\n", &n);
	printf("read n= %d\n", n);

	if (n != N) {
		printf("Restart file is not compatible with current system.\n");
		exit(1);
	}

	for (i=0; i<2*n; i++)
		fscanf (restart, "%f\t%f\t%f\n", &x[i][0],&x[i][1],&x[i][2]);
	
	fscanf(restart, "%c\n", &check);

	if (check != 'v') {
		printf("Something wrong with restart file.\n");
		exit(1);
	}

	for (i=0; i<2*n; i++)
		fscanf (restart, "%f\t%f\t%f\n", &v[i][0],&v[i][1],&v[i][2]);
	
	fscanf(restart, "%c\n", &check);

	if (check != 'f') {
		printf("Something wrong with restart file.\n");
		exit(1);
	}



	for (i=0; i<2*n; i++)
		fscanf (restart, "%f\t%f\t%f\n", &f[i][0],&f[i][1],&f[i][2]);

	fclose (restart);
}

static void writeRestart(char * filename) {
	int i;

	FILE *restart = fopen(filename, "w");

	fprintf (restart, "%d\n", N);

	for (i=0; i<2*N; i++)
		fprintf (restart, "%f\t%f\t%f\n", x[i][0],x[i][1],x[i][2]);
	
	fprintf(restart, "v\n");

	for (i=0; i<2*N; i++)
		fprintf (restart, "%f\t%f\t%f\n", v[i][0],v[i][1],v[i][2]);
	
	fprintf(restart, "f\n");

	for (i=0; i<2*N; i++)
		fprintf (restart, "%f\t%f\t%f\n", f[i][0],f[i][1],f[i][2]);

	fclose (restart);
}


