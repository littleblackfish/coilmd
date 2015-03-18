#! /usr/bin/env python2

from MDAnalysis import *
from numpy import *
from sys import argv
#from pylab import *

def dot(a,b) : return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

print 'Loading trajectory, might take a while.'
system=Universe(argv[1], argv[2])

P1= system.selectAtoms('resname dna1') 	#select atoms from 1 strand

maxtau=40 				#cutoff for tau

n=len(P1)				#final length of chain

nframes=system.trajectory.numframes-1	#length of trajectory

tautop,taubottom,correl=zeros(maxtau),zeros(maxtau),zeros(maxtau)

print 'Calculating correlations..'   

for y in range(nframes) :
	print  "%.1f percent complete  \r " %(y*100./nframes) ,
	
	tmp=P1.coordinates()
	r=zeros([(n-1) , 3])
	for x in range(n-1) : 		#building r vectors along backbone	
		r[x]= tmp[x+1]-tmp[x] 

	for tau in range(maxtau):	#calculating dot products
		for x in range (n-tau-1):
			tautop[tau]+=dot(r[x+tau],r[x])
			taubottom[tau]+=dot(r[x],r[x])

	system.trajectory.next()

for x in range(len(tautop)) :		#calculating correlations
	correl[x] = tautop[x]/taubottom[x]

s=array([arange(maxtau), correl])

savetxt('persistence', s.T)  
