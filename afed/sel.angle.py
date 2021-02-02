#!/usr/bin/env python3

import time
import numpy as np
import argparse

parser = argparse.ArgumentParser(
    description="*** Train a model. ***")
parser.add_argument('-t','--threshold', default = 5, type=float,
                    help='the input angles ')
parser.add_argument('INPUT', help='the input angles ')
parser.add_argument('OUTPUT', help='the output angles ')
args = parser.parse_args()

data = np.loadtxt (args.INPUT)
sel  = np.array([data[0]])
sel_idx = np.array([0], dtype=np.int)
threshold = args.threshold

nframes=data.shape[0]
nvalues=data.shape[1]

print (sel)

def diff_pbc (a, b) :
    nvalue = len(a)
    diff = a - b
    for ii in range(nvalue) :
        if diff[ii] >= 180 : 
            diff[ii] -= 360
        elif diff[ii] < -180 :
            diff[ii] += 360
    return diff

for ii in range (1, nframes) :
    do_append = True
    sel_range = len(sel)
    for jj in range(sel_range) :
        diff = diff_pbc (data[ii], sel[jj])
        if np.linalg.norm(diff) < threshold : 
            do_append = False
            break
    if do_append : 
        sel = np.reshape (np.append (sel, data[ii]), [-1,nvalues])
        sel_idx = np.append (sel_idx, ii)
        print ("selected: %d centers" % sel_range)
        # print (sel)

np.savetxt (args.OUTPUT, sel, fmt = "%.6f")
np.savetxt (args.OUTPUT+".idx", sel_idx, fmt = "%d")
