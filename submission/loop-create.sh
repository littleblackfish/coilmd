#!/bin/bash 

# bash script to generate a new branch in the folder hierarchy
# input is N
# run from data root.

N=$1

mkdir $N;
cd $N;

mkdir {cir,lin};

for type in ???; do
	cd $type;
	mkdir 0.{10..30};
	cd ..;
done


