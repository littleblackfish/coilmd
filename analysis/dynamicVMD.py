#!/usr/bin/env python2

from sys import argv

n=int(argv[1])

print 'mol new traj.vtf type vtf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all'
print 'mol delrep 0 top'
print 'mol representation DynamicBonds 1.2 0.05 6.'
print 'mol color colorID 6'
print 'mol material Transparent'

for  x in range(0,2*n,2) :
    print 'mol selection {index ' +repr(x) + ' ' +repr(x+1) + '}'
    print 'mol addrep top'

print 'mol representation Licorice 0.100000 10.000000 10.000000'
print 'mol material Opaque'
print 'mol color Chain'
print 'mol selection {chain 0}'
print 'mol addrep top'
print 'mol selection {chain 1}'
print 'mol addrep top'



