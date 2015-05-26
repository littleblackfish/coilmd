#!/bin/bash

# bulk submission script 
# input is sge submission file to be used and N

xsubmitfile=$1

for N in $2; do 
	cd $N; 
	for type in ???; do 
		cd $type; 
		for temp in ?.??; do 
			cd $temp; 
			qsub -N ${type}${N}-${temp} ${xsubmitfile} $type $N ${temp}; 
			cd ..; 
		done ; 
		cd ..; 
	done; 
	cd .. ; 
done
