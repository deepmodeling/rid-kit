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
    with tf.gfile.GFile(frozen_graph_filename, "rb") as f:
        graph_def = tf.GraphDef()
        graph_def.ParseFromString(f.read())
    with tf.Graph().as_default() as graph:
        tf.import_graph_def(
            graph_def, 
            input_map=None, 
            return_elements=None, 
            name=prefix, 
            producer_op_list=None
        )
    return graph

def test_ef (sess, xx) :
    graph = sess.graph

    inputs  = graph.get_tensor_by_name ('load/inputs:0')
    o_energy= graph.get_tensor_by_name ('load/o_energy:0')
    o_forces= graph.get_tensor_by_name ('load/o_forces:0')
    
    zero4 = np.zeros (xx.shape)
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
    return np.array (stds)

def compute_model_dev (models, data) :
    nframes = data.shape[0]
    cv_dim = data.shape[1]
    forces = []
    for ii in models :
        graph = load_graph (ii)
        with tf.Session(graph = graph) as sess:        
            ee, ff = test_ef (sess, data)
            forces = np.append (forces, ff)
    forces = np.reshape (forces, [len(models), nframes, cv_dim])
    forces *= f_cvt

    # std_f = np.std(forces, axis = 0)
    # std_f = np.linalg.norm(std_f, axis = 1)
    std_f = compute_std (forces)

    assert (std_f.shape[0] == nframes)

    return std_f
