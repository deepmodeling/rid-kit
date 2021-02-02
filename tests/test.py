#!/usr/bin/env python3

import re
import os
import sys
import argparse
import numpy as np
import tensorflow as tf

kbT = (8.617343E-5) * 300 
beta = 1.0 / kbT
f_cvt = 96.485

def load_graph(frozen_graph_filename, 
               prefix = 'load'):
    # We load the protobuf file from the disk and parse it to retrieve the 
    # unserialized graph_def
    with tf.gfile.GFile(frozen_graph_filename, "rb") as f:
        graph_def = tf.GraphDef()
        graph_def.ParseFromString(f.read())

    # Then, we can use again a convenient built-in function to import a graph_def into the 
    # current default Graph
    with tf.Graph().as_default() as graph:
        tf.import_graph_def(
            graph_def, 
            input_map=None, 
            return_elements=None, 
            name=prefix, 
            op_dict=None, 
            producer_op_list=None
        )
    return graph

def test_ef (sess, xx) :
    graph = sess.graph

    inputs  = graph.get_tensor_by_name ('load/inputs:0')
    o_energy= graph.get_tensor_by_name ('load/o_energy:0')
    o_forces= graph.get_tensor_by_name ('load/o_forces:0')

    cv_dim = xx.shape[1]
    zero4 = np.zeros ([xx.shape[0], cv_dim])
    data_inputs = np.concatenate ((xx, zero4), axis = 1)
    feed_dict_test = {inputs: data_inputs}

    data_ret = sess.run ([o_energy, o_forces], feed_dict = feed_dict_test)
    return data_ret[0], data_ret[1]


def compute_std (forces) :
    nmodels = forces.shape[0]
    nframes = forces.shape[1]
    ncomps = forces.shape[2]
    
    stds = []
    for ii in range (nframes) :
        # print ( forces[0, ii], forces[1, ii], forces[2, ii], forces[3, ii])
        avg_std = 0
        for jj in range (ncomps) :
            mystd = np.std (forces[:, ii, jj])
            avg_std += mystd * mystd
        avg_std = np.sqrt (avg_std / float(ncomps))
        stds.append (avg_std)
    
    # avg_f = np.average (forces, axis = 0)
    # norm_force = np.linalg.norm(avg_f, axis = 1)
    # print (norm_force)
    # print (np.average (norm_force))
    # print (stds)
    # print (np.max(stds), np.average(stds), np.std(stds))
    # rel = stds / norm_force
    # print (rel)
    # print (np.max(rel), np.average(rel), np.std(rel))
    return np.array (stds)
    

def _main () :
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--models", default=[], nargs = '*', type=str, 
                        help="Frozen models file to test")
    parser.add_argument("-d", "--data", type=str, 
                        help="The data for test")
    parser.add_argument("-D", "--cv-dim", type=int,
                        help="The dim of CV")
    args = parser.parse_args()

    models = args.models    
    cv_dim = args.cv_dim
    data_ = np.loadtxt (args.data)
    assert (data_.shape[1] >= cv_dim), "the data size should be larger than the CV dim size"
    data = data_[:,:cv_dim]
    nframes = data.shape[0]
    
    energies = []
    forces = []
    for ii in models :
        graph = load_graph (ii)
        with tf.Session(graph = graph) as sess:        
            ee, ff = test_ef (sess, data)
            energies = np.append (energies, ee)
            forces = np.append (forces, ff)
            # if len(forces) == 0 :
            #     forces = np.array ([ff])
            # else :
            #     forces = np.concatenate (forces, np.array([ff]))

    energies = np.reshape (energies, [len(models), nframes, 1])
    forces = np.reshape (forces, [len(models), nframes, cv_dim])
    energies *= f_cvt
    forces *= f_cvt
    for ii in range(len(models)) :        
        energies[ii] -= np.max (energies[ii])
    norm_f = np.linalg.norm (forces, axis = 2)
    
    avg_e_std = compute_std (energies)
    avg_f_std = compute_std (forces)
    rel_f_std = np.divide (avg_f_std, norm_f[0])

    '''
    print (energies)
    ae = np.average((energies), axis = 0)
    ad = np.std((energies), axis = 0) / np.sqrt(len(models))
    res = np.concatenate ((ae - ae[0], ad), axis = 1)
    print (res)
    '''

    print ("max_e:     %.2e  min_e:    %.2e  diff:      %.2e" % (np.max(energies), np.min(energies), np.max(energies) - np.min(energies)))
    print ("max_|f|:   %.2e  min_|f|:   %.2e  avg_|f|:   %.2e" % (np.max(norm_f), np.min(norm_f), np.average(norm_f)))
    print ("max_std_e: %.2e  min_std_e: %.2e  avg_std_e: %.2e" % (np.max(avg_e_std), np.min(avg_e_std), np.average(avg_e_std)))
    # print ("max_f    : %.2e  min f    :%.2e  diff     : %.2e" % (np.max(forces), np.min(forces), np.max(forces) - np.min(forces)))
    print ("max_std_f: %.2e  min_std_f: %.2e  avg_std_f: %.2e" % (np.max(avg_f_std), np.min(avg_f_std), np.average(avg_f_std)))    
    print ("max_rel_f: %.2e  min_rel_f: %.2e  avg_rel_f: %.2e" % (np.max(rel_f_std), np.min(rel_f_std), np.average(rel_f_std)))    

if __name__ == '__main__':
    _main()
