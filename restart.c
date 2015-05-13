int readRestart(char * filename) {
	int i,n;
	char check;
	FILE *restart = fopen(filename, "r");

	if (restart) {
		
		printf("Restart file found.\n");
		
#ifdef _OPENMP
		fscanf(restart,"%d\t%u\n", &n, &seed[0]);
#else
		fscanf(restart,"%d\t%u\n", &n, &seed);
#endif

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
		return 1;
	}
	else { 
		printf("Restart file NOT found.\n");
		return 0;
	}
}

static void writeRestart(char * filename) {
	int i;

	FILE *restart = fopen(filename, "w");

#ifdef _OPENMP
	fprintf (restart, "%d\t%u\n", N, seed[0]);
#else
	fprintf (restart, "%d\t%u\n", N, seed);
#endif

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
