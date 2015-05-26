#! /usr/bin/env python2

# script to calculate the linking number of a circular system

from sys import argv,stdout
from MDAnalysis import *
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

# plots each 3d curves in a list with a different color 

def plot3d (R) :
	fig=figure(1)
	ax=fig.gca(projection='3d')
	
        for r in R :
            ax.plot(r.T[0],r.T[1],r.T[2])

	show()

# test case : two interlocked circles should have link = 1

def circles (npoints) :
    theta= linspace(0, 2*pi, npoints+1) 
    theta = theta[0:-1]
    circle1 =array( [cos(theta), sin(theta), zeros(npoints)]).T
    circle2 =array( [zeros(npoints), 1+sin(theta), cos(theta)]).T
    
    print link(circle1, circle2)

    plot3d([circle1, circle2])


import pyximport ; pyximport.install()
from hasan import link
from hasan import resample


if __name__ == '__main__': 
    
    system=Universe(argv[1], argv[2])
    
    P1= system.selectAtoms('resname dna1') 	#select atoms from 1st strand
    P2= system.selectAtoms('resname dna2') 	#select atoms from 2nd strand
    nframes=system.trajectory.numframes-1	#length of trajectory
    n=len(P1)				        #length of chain

    lnlist=zeros(nframes+1)


    for t in [nframes] :  #in range(nframes+1) :
#    for t in range(nframes+1) :  #in range(nframes+1) :
#    for t in range(85280,85310) :  #in range(nframes+1) :
#    for t in range(101, 110):
        system.trajectory[t]

    	r1 = P1.coordinates()
        r2 = P2.coordinates()

        r1=resample(r1,3*n)
        r2=resample(r2,3*n)
        
        tmp = link(r1,r2)
       # print t,tmp
        print tmp

        lnlist[t] = tmp
    
    
    
    savetxt('ln.dat', lnlist)  
