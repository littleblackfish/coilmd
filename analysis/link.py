#! /usr/bin/env python2

# script to calculate the linking number of a circular system

from sys import argv
from MDAnalysis import *
from numpy import *
from scipy import interpolate
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

# some local functions to reduce overhead

def dot(a,b) : return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def norm(a) : return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])

# plots each 3d curves in a list with a different color 

def plot3d (R) :
	fig=figure(1)
	ax=fig.gca(projection='3d')
	
        for r in R :
            ax.plot(r.T[0],r.T[1],r.T[2])

	show()

# interpolates a curve r with n points in total

def resample(r, n) :
        r = vstack( [r,r[0]] )
        tck, u = interpolate.splprep(r.T, s=0)
        return array(interpolate.splev(linspace(0,1,n+1),tck)).T[:-1]

# calculates derivative of a given curve at each point

def derivative (r) :
	n=len(r)
	dr=empty((n,3))
	
	for i in range(n) :
		dr[i]= 0.5* ( r[(i+1)%n] - r[i-1] )
	
	return dr

# calculates the linking number of given closed curve

def link(r1, r2) :
    ln = 0
    n  = len(r1)

    dr1 = derivative(r1)
    dr2 = derivative(r2)
    
    for i in range(n)  :
        for j in range(n) :
            r12 = r1[i]-r2[j]
            ln += dot ( r12 , cross(dr1[i], dr2[j]) ) / (norm(r12)**3)


    return ln /(4*pi)

# test case : two interlocked circles should have link = 1

def circles (npoints) :
    theta= linspace(0, 2*pi, npoints+1) 
    theta = theta[0:-1]
    circle1 =array( [cos(theta), sin(theta), zeros(npoints)]).T
    circle2 =array( [zeros(npoints), 1+sin(theta), cos(theta)]).T
    
    print link(circle1, circle2)

    plot3d([circle1, circle2])


if __name__ == '__main__': 
    
    print 'Loading trajectory, might take a while.'
    system=Universe(argv[1], argv[2])
    
    P1= system.selectAtoms('resname dna1') 	#select atoms from 1st strand
    P2= system.selectAtoms('resname dna2') 	#select atoms from 2nd strand
    
    
    nframes=system.trajectory.numframes-1	#length of trajectory
    
    print 'Calculating linking number..'   
    n=len(P1)				#final length of chain
    
    lnlist=zeros(nframes)
    
    for t in range(nframes) :
    	print  "%.1f percent complete  \r " %(t*100./nframes) ,
    	
    	r1 = P1.coordinates()
        r2 = P2.coordinates()

        r1=resample(r1,500)
        r2=resample(r2,500)

        lnlist[t] = link(r1,r2)
    
    	system.trajectory.next()
    
    savetxt('ln.dat', lnlist)  
