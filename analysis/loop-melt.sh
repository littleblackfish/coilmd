#!/bin/bash

# This is a loop that calculates a melting curve for a given N and type 

N=$1
type=$2

cd $N/${type}
rm melt.dat

coilmd="~/coilmd"
command="${coilmd}/analysis/melt.awk "

for temp in ?.??; do 
	cd ${temp}; 
	pwd; 
	awk 'NR%10 == 0' bubbles.dat > bubbles-sparse.dat; 
	rate=$(~/coilmd/analysis/melt.awk bubbles-sparse.dat);
	echo ${temp} ${rate} >> ../melt.dat ; 
	cd .. ; 
done
