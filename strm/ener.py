#!/usr/bin/env python3

import os
import argparse
import numpy as np
import subprocess as sp
import libstring as libs

def main () :
    parser = argparse.ArgumentParser(
        description="*** Compute the free energy along a string. ***")

    parser.add_argument('string_folder', type = str,
                        help='The folder of the string to be computed.')
    parser.add_argument('-o', '--output', default = "energy.out",
                        help='The output energy in string folder.')

    args = parser.parse_args()

    string_file = args.string_folder + "/string.out"
    string_force = args.string_folder + "/force.out"
    if not os.path.exists (string_file) :
        raise RuntimeError ("cannot find string file " + string_file)
    if not os.path.exists (string_force) :
        raise RuntimeError ("cannot find string force file " +
                            string_force +
                            ". maybe the simulation is not finished")

    string = np.loadtxt (string_file)
    force = np.loadtxt (string_force)
    numb_node = string.shape[0]
    dim = string.shape[1]    

    # compute arc (alpha)
    alpha_seg = libs.arc_seg (string)
    alpha = libs.arc_norm (string)

    # integrate the energy
    energy = np.zeros (alpha.shape)
    tagent = libs.compute_string_tegent (alpha, string)
    for ii in range (1, numb_node) :
        v0 = np.dot (tagent[ii-1], force[ii-1])
        v1 = np.dot (tagent[ii]  , force[ii]  )
        vn = force[ii] - v1 * tagent[ii]
        # print (str(v1) +
        #        "  t: " + str(v1 * tagent[ii]) +
        #        "  |t|:  " + str(np.sqrt(np.dot(v1 * tagent[ii], v1 * tagent[ii]))) +
        #        "  n: " + str(vn) +
        #        "  |n|: " + str(np.sqrt(np.dot(vn,vn))) +
        #        "  t.n: " + str(np.dot(v1 * tagent[ii], vn))
        #        )
        de = 0.5 * (v0 + v1) * alpha_seg[ii]
        energy[ii] = energy[ii-1] - de

    new_shape = np.array ([alpha.shape[0], 2])
    result = np.zeros (new_shape)
    result[:,0] = alpha
    result[:,1] = energy
    string_energy = args.string_folder + "/" + args.output
    np.savetxt (string_energy, result)
    
if __name__ == "__main__":
    main ()
    
