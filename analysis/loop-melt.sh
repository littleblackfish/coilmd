#!/bin/bash

# This is a loop that calculates a melting curve for a given N and type 

N=$1
type=$2

cd $N/${type}

$coilmd="~/scratch/coilmd"

for x in ?.??; do 
	cd $x; 
	pwd; 
	awk 'NR%10 == 0' bubbles.dat > bubbles-sparse.dat; 
	echo $x $(${coilmd}/analysis/melt.awk bubbles-sparse.dat) >> ../melt.dat ; 
	cd .. ; 
done
