#!/bin/bash

# bulk submission script 

$xsubmit="../../../xsubmit-somon "


for N in *; do 
	cd $N; 
	for type in ???; do 
		cd $type; 
		for temp in ?.??; do 
			cd $temp; 
			qsub -N ${type}${N}-$temp $xsubmit 
			$type $N ${temp}; 
			cd ..; 
		done ; 
		cd ..; 
	done; 
	cd .. ; 
done
