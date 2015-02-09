#!/bin/awk -f

# awk script to calculate the dissociation ratio quickly

{
	s=0;
	for (i=1; i<=NF; i++) 	s+=$i
	
	molten[NR] = 1-(s/NF)
}

END{
	# binning
	
	binCount=10;
	binSize = int(NR/binCount);
	currentBin=1;

	for (i=1; i <= NR; i++ ) {
		if (i%binSize ==0) { 
			bin[currentBin] /= count;
			currentBin ++;
			count = 0;
		}
		bin[currentBin] += molten[i];
		count ++;
	}
	
	# mean over bins
	
	for (i=1; i<=binCount; i++)  mean += bin[i];

	mean /= binCount;

	# std over bins

	for (i=1; i<=binCount; i++) std += (bin[i]-mean)**2;

	std = sqrt(std/(binCount-1))
		
	err = std/sqrt(binCount)

	print mean, err
}

