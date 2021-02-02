#!/usr/bin/env python3

import numpy as np; 

kk=np.loadtxt('kappa.out'); 
cc=np.loadtxt('centers.out'); 
data=np.loadtxt ('avgins.centers.out'); 

avgins=data[:,0]; 
erravg=data[:,1];

diff=np.zeros (avgins.shape)

for ii in range(len(avgins)) :
    diff[ii] = avgins[ii] - cc[ii]
    if   diff[ii] >= np.pi :
        diff[ii] -= np.pi * 2.
    elif diff[ii] < -np.pi :
        diff[ii] += np.pi * 2.

ff = np.multiply ( kk, diff ); 
fe = np.multiply ( kk, erravg); 
np.savetxt ('force.out',  ff, fmt='%.10e');
np.savetxt ('ferror.out', fe, fmt='%.10e');
