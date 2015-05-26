from numpy import *
from scipy import interpolate

cdef float  dot(float a[3], float b[3]) : return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
cdef float  norm(float a[3]) : return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])

def derivative (r) :
    n=len(r)
    dr=empty((n,3), dtype=float32)
    cdef int i
    for i in range(n) :
        dr[i]= 0.5* ( r[(i+1)%n] - r[i-1] )
    return dr

def link(r1, r2) :
    dr1 = derivative(r1)
    dr2 = derivative(r2)

    cdef int i,j,n
    cdef float r12[3]
    cdef float drcross[3]
    cdef float ln=0
    
    n  = len(r1)
    
    for i in range(n)  :
        for j in range(n) :
            r12 = r1[i]-r2[j]
            drcross = cross(dr1[i],dr2[j])
            ln += dot ( r12 , drcross ) / (norm(r12)**3)
    
    return ln /(4*pi)


def resample(r, n) :
    r = vstack( [r,r[0]] )
    tck, u = interpolate.splprep(r.T, s=0)
    return array(interpolate.splev(linspace(0,1,n+1),tck), dtype=float32).T[:-1]
