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

def graph_op_name (graph) :
    names = []
    for op in graph.get_operations():
        print (op.name)

def test_ef (sess, xx) :
    graph = sess.graph

    inputs  = graph.get_tensor_by_name ('load/inputs:0')
    o_energy= graph.get_tensor_by_name ('load/o_energy:0')
    o_forces= graph.get_tensor_by_name ('load/o_forces:0')

    zero4 = np.zeros ([xx.shape[0], 4])
    data_inputs = np.concatenate ((xx, zero4), axis = 1)
    feed_dict_test = {inputs: data_inputs}

    data_ret = sess.run ([o_energy, o_forces], feed_dict = feed_dict_test)
    return data_ret[0], data_ret[1]

def test_e (sess, xx) :
    graph = sess.graph

    inputs  = graph.get_tensor_by_name ('load/inputs:0')
    o_energy= graph.get_tensor_by_name ('load/o_energy:0')

    zero4 = np.zeros ([xx.shape[0], 4])
    data_inputs = np.concatenate ((xx, zero4), axis = 1)
    feed_dict_test = {inputs: data_inputs}

    data_ret = sess.run ([o_energy], feed_dict = feed_dict_test)
    return data_ret[0]

def value_array_pp (sess, ngrid) :
    delta  = 2.0 * np.pi / ngrid
    delta2 = delta * delta
    xx = np.linspace(-np.pi + 0.5*delta, np.pi - 0.5 * delta, (ngrid))
    yy = np.linspace(-np.pi + 0.5*delta, np.pi - 0.5 * delta, (ngrid))
    zz = np.zeros([len(xx) * len(yy)])

    my_grid    = np.zeros((ngrid*ngrid,2))
    zero_grid  = np.zeros((ngrid*ngrid,1))
    for i in range(ngrid):
        for j in range(ngrid):
            my_grid[i*ngrid + j, 0] = xx[i]
            my_grid[i*ngrid + j, 1] = yy[j]

    for i in range(len(xx)):
        print("computing grid: %d" % i)
        for j in range(len(yy)):
            ve = test_e (sess, np.concatenate((zero_grid + xx[i], zero_grid + yy[j], my_grid), axis = 1))
            zz[i*(ngrid)+j] = (kbT) * np.log(np.sum(delta2 * np.exp(-beta * ve)))

    zz = zz - np.max(zz)
    return xx, yy, zz

def value_array_pt (sess, ngrid) :
    delta  = 2.0 * np.pi / ngrid
    delta2 = delta * delta
    xx = np.linspace(-np.pi + 0.5*delta, np.pi - 0.5 * delta, (ngrid))
    yy = np.linspace(-np.pi + 0.5*delta, np.pi - 0.5 * delta, (ngrid))
    zz = np.zeros([len(xx) * len(yy)])

    my_grid    = np.zeros((ngrid*ngrid,2))
    zero_grid  = np.zeros((ngrid*ngrid,1))
    for i in range(ngrid):
        for j in range(ngrid):
            my_grid[i*ngrid + j, 0] = xx[i]
            my_grid[i*ngrid + j, 1] = yy[j]

    for i in range(len(xx)):
        print("computing grid: %d" % i)
        for j in range(len(yy)):
            ve = test_e (sess, np.concatenate((zero_grid + xx[i], np.reshape(my_grid[:,0], [-1,1]), zero_grid + yy[j], np.reshape(my_grid[:,1], [-1,1])), axis = 1))
            zz[i*(ngrid)+j] = (kbT) * np.log(np.sum(delta2 * np.exp(-beta * ve)))

    zz = zz - np.max(zz)
    return xx, yy, zz

def value_array_tz (sess, ngrid) :
    delta  = 2.0 * np.pi / ngrid
    delta2 = delta * delta
    xx = np.linspace(-np.pi + 0.5*delta, np.pi - 0.5 * delta, (ngrid))
    yy = np.linspace(-np.pi + 0.5*delta, np.pi - 0.5 * delta, (ngrid))
    zz = np.zeros([len(xx) * len(yy)])

    my_grid    = np.zeros((ngrid*ngrid,2))
    zero_grid  = np.zeros((ngrid*ngrid,1))
    for i in range(ngrid):
        for j in range(ngrid):
            my_grid[i*ngrid + j, 0] = xx[i]
            my_grid[i*ngrid + j, 1] = yy[j]

    for i in range(len(xx)):
        print("computing grid: %d" % i)
        for j in range(len(yy)):
            ve = test_e (sess, np.concatenate((my_grid, zero_grid + xx[i], zero_grid + yy[j]), axis = 1))
            zz[i*(ngrid)+j] = (kbT) * np.log(np.sum(delta2 * np.exp(-beta * ve)))

    zz = zz - np.max(zz)
    return xx, yy, zz

def value_array_pz (sess, ngrid) :
    delta  = 2.0 * np.pi / ngrid
    delta2 = delta * delta
    xx = np.linspace(-np.pi + 0.5*delta, np.pi - 0.5 * delta, (ngrid))
    yy = np.linspace(-np.pi + 0.5*delta, np.pi - 0.5 * delta, (ngrid))
    zz = np.zeros([len(xx) * len(yy)])

    my_grid    = np.zeros((ngrid*ngrid,2))
    zero_grid  = np.zeros((ngrid*ngrid,1))
    for i in range(ngrid):
        for j in range(ngrid):
            my_grid[i*ngrid + j, 0] = xx[i]
            my_grid[i*ngrid + j, 1] = yy[j]

    for i in range(len(xx)):
        print("computing grid: %d" % i)
        for j in range(len(yy)):
            ve = test_e (sess, np.concatenate((np.reshape(my_grid[:,0], [-1,1]), zero_grid + xx[i], np.reshape(my_grid[:,1], [-1,1]), zero_grid + yy[j]), axis = 1))
            zz[i*(ngrid)+j] = (kbT) * np.log(np.sum(delta2 * np.exp(-beta * ve)))

    zz = zz - np.max(zz)
    return xx, yy, zz

def print_array (fname, xx, yy, zz0, zz1, zz2, zz3) :
    with open(fname, 'w') as fp :
        lx = len(xx)
        for ii in range (len(xx)) :
            for jj in range (len(yy)) :
                fp.write ("%f %f %f %f %f %f\n" % (xx[ii], yy[jj], zz0[ii*lx+jj], zz1[ii*lx+jj], zz2[ii*lx+jj], zz3[ii*lx+jj]) )
            fp.write ("\n")

def _main () :
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--model", default=["frozen_model.pb"], type=str, nargs="*",
                        help="Frozen model file to import")
    parser.add_argument("-n", "--numb-grid", default=60, type=int, 
                        help="The number of data for test")
    parser.add_argument("-o", "--output", default="fe.out", type=str, 
                        help="output free energy")
    args = parser.parse_args()
    
    count = 0
    for ii in args.model :
        graph = load_graph(ii)
        with tf.Session(graph = graph) as sess:        
            xx, yy, zz0 = value_array_pp (sess, args.numb_grid)
            xx, yy, zz1 = value_array_tz (sess, args.numb_grid)
            xx, yy, zz2 = value_array_pt (sess, args.numb_grid)
            xx, yy, zz3 = value_array_pz (sess, args.numb_grid)
            zz0 *= f_cvt
            zz1 *= f_cvt
            zz2 *= f_cvt
            zz3 *= f_cvt
            if count == 0: 
                avg0 = zz0
                avg1 = zz1
                avg2 = zz2
                avg3 = zz3
            else :
                avg0 += zz0
                avg1 += zz1
                avg2 += zz2
                avg3 += zz3
        count += 1
    avg0 /= float(count)
    avg1 /= float(count)
    avg2 /= float(count)
    avg3 /= float(count)

    avg0 -= np.max(avg0)
    avg1 -= np.max(avg1)
    avg2 -= np.max(avg2)
    avg3 -= np.max(avg3)
    print_array(args.output, xx, yy, avg0, avg1, avg2, avg3)

if __name__ == '__main__':
    _main()
