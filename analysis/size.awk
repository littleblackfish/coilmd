#!/bin/awk -f

# awk script to calculate the dissociation ratio quickly

{
	inside = 0;
       	size = 0;

	for (i=1; i<=NF; i++) 
	{
		if (inside) {
			if ($i == 0) size ++ ;
			else {
				sizeCount[size] ++;
				inside = 0;
				size = 0;
			}
		}
		else {
			if ($i ==1) {
				inside = 1;
				size ++;
			}
		}

	} 
	
	}; 

END{
	for (i in sizeCount )
		print i, sizeCount[i]
}

