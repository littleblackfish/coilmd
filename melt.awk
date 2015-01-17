#!/bin/awk -f

# awk script to calculate the dissociation ratio quickly

{for (i=1; i<=NF; i++) s=s+$i}; 

END{print 1-s/NR/NF}

