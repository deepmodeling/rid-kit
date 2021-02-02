#!/usr/bin/env python3

import time
import numpy as np
import argparse

parser = argparse.ArgumentParser(
    description="*** Train a model. ***")
parser.add_argument('INPUT', help='the input plumed output ')
parser.add_argument('OUTPUT', help='the input plumed output ')
args = parser.parse_args()

data = np.loadtxt (args.INPUT)

nframes=data.shape[0]
nvalues=data.shape[1]

for ii in range (1, nframes) :
    for jj in range (1, nvalues) :
        if   data[ii,jj] - data[0,jj] >= np.pi :
            data[ii,jj] -= np.pi * 2.
        elif data[ii,jj] - data[0,jj] < -np.pi :
            data[ii,jj] += np.pi * 2.

np.savetxt (args.OUTPUT, data, fmt="%.6f")
