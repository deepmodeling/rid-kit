#!/usr/bin/env python3

import re
import os
import sys
import argparse
import logging
import numpy as np
from typing import Union, List
try:
    import tensorflow.compat.v1 as tf
    tf.disable_v2_behavior()
except ImportError:
    import tensorflow as tf
from rid.common.tensorflow.graph import load_graph
from rid.constants import kbT, beta, f_cvt


logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def test_ef(sess, data_in):
    graph = sess.graph

    inputs = graph.get_tensor_by_name('load/inputs:0')
    o_energy = graph.get_tensor_by_name('load/o_energy:0')
    o_forces = graph.get_tensor_by_name('load/o_forces:0')
    drop_out_rate = graph.get_tensor_by_name('load/drop_out_rate:0')

    zero4 = np.zeros([data_in.shape[0], data_in.shape[1]])
    data_inputs = np.concatenate((data_in, zero4), axis=1)
    feed_dict_test = {inputs: data_inputs, drop_out_rate: 0.0}

    data_ret = sess.run([o_energy, o_forces], feed_dict=feed_dict_test)
    return data_ret[0], data_ret[1]


def compute_std(forces):
    stds = np.mean( np.std(forces, axis=0) ** 2, axis=-1 ) ** 0.5
    return stds


def make_std(
        data: Union[List, np.ndarray], 
        models: List = ["graph.000.pb"]
        # threshold: float = 1.0
    ):

    nframes = data.shape[0]

    forces = []
    for model in models:
        graph = load_graph(str(model))
        with tf.Session(graph=graph) as sess:
            _, force = test_ef(sess, data)
            forces = np.append(forces, force)

    forces = np.reshape(forces, [len(models), nframes, -1])
    forces *= f_cvt

    avg_std = compute_std(forces)
    return avg_std

