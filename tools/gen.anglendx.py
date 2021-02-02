#!/usr/bin/env python3

import argparse
import numpy as np

# numb_calpha = 10
# angle_atoms = [ [5,7,9,15],
#                 [7,9,15,17],
#                 [6,5,7,9],
#                 [9,15,17,18] ]

# def make_angle_def (angle_name,
#                     angle_atoms,
#                     idx_cum,
#                     alpha_idx) :
#     mylist = ""
#     for ii in range (0,len(angle_atoms)) :
#         mylist += " " + str(idx_cum[alpha_idx] + angle_atoms[ii])
#     shifted_name = (prt_shift_fmt % alpha_idx) + "-" + angle_name
#     return ("[ " + shifted_name + " ]\n" +
#             mylist)

def make_angle_def (idx_cum,
                    angle_idxes,
                    fmt_alpha,
                    fmt_angle) :
    ret = ""
    for ii in range(len(idx_cum)) :
        for jj in range(len(angle_idxes)) :
            angle_print = (fmt_alpha % ii) + "-" + (fmt_angle % jj)
            mylist = ""
            for kk in angle_idxes[jj] :
                mylist += " " + str(idx_cum[ii] + kk)
            ret += "[ " + angle_print + " ]\n" 
            ret += mylist
            ret += "\n"
    return ret

def make_shift_cum (numb_calpha) :
    idx_shift = np.ones(numb_calpha, dtype=np.int) * 10
    idx_shift[0] = 0
    idx_cum = np.cumsum (idx_shift)
    return idx_cum

def _main () :
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--numb-alpha", default=0, type=int, 
                        help="the number of alpha carbons")
    parser.add_argument("-i", "--angle-indexes", default=[], nargs= '*', type=int, 
                        help="the indexes of dihedral angle atoms")
    parser.add_argument("--alpha-fmt", default="%02d", type=str, 
                        help="format of alpha")
    parser.add_argument("--angle-fmt", default="%02d", type=str, 
                        help="format of angle")
    
    args = parser.parse_args()
    
    numb_calpha = args.numb_alpha
    # shape should be numb_angles x 4
    angle_idxes = np.reshape (np.array(args.angle_indexes, dtype=np.int), [-1,4])    

    if numb_calpha == 0: 
        return

    idx_cum = make_shift_cum (numb_calpha)
    
    print (make_angle_def(idx_cum, angle_idxes, args.alpha_fmt, args.angle_fmt))

if __name__ == '__main__':
    _main()

