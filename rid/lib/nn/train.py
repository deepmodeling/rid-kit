#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import tensorflow as tf
import numpy as np
from rid.lib.nn.model import Reader, Model


class Config(object):
    def __init__(self, cv_dim_list):
        self.batch_size = 64
        self.n_displayepoch = 200
        # use Batch Normalization or not
        self.useBN = False  # True
        self.num_epoch = 3000
        self.use_mix = False
        self.old_ratio = 7
        self.n_neuron = [240, 120, 60, 30]
        self.starter_learning_rate = 0.003
        self.decay_steps = 10
        self.decay_steps_inner = 0
        self.decay_rate = 0.96
        self.data_path = './data/'
        self.restart = False
        self.resnet = False
        self.graph_file = None
        self.cv_dim_list = cv_dim_list
        self.cv_dim = sum(self.cv_dim_list)
        self.display_in_training = True


def reset_batch_size(config):
    tr_data = np.loadtxt(config.data_path+'data.raw')
    if tr_data.shape[0] < config.batch_size:
        config.batch_size = tr_data.shape[0]
        print("using new batch_size of %d" % config.batch_size)


def print_conf(config, nthreads):
    print("# restart           " + str(config.restart))
    print("# num_threads       %d" % nthreads)
    print("# neurons           " + str(config.n_neuron))
    print("# batch size        " + str(config.batch_size))
    print("# num_epoch         " + str(config.num_epoch))
    print("# use_mix           " + str(config.use_mix))
    print("# old_ratio         " + str(config.old_ratio))
    print("# lr_0              " + str(config.starter_learning_rate))
    print("# decay_steps       " + str(config.decay_steps))
    print("# decay_steps_inner " + str(config.decay_steps_inner))
    print("# decay_rate        " + str(config.decay_rate))
    print("# resnet            " + str(config.resnet))
    print("# graph_file        " + str(config.graph_file))


def set_conf(cv_dim_list,
             neurons=[240, 120, 60, 30],
             numb_threads=8,
             resnet=True,
             use_mix=False,
             restart=False,
             batch_size=128,
             epoches=12000,
             lr=0.0008,
             decay_steps=120,
             decay_rate=0.96,
             old_ratio=7.0,
             decay_steps_inner=0):
    config = Config(cv_dim_list)
    config.n_neuron = neurons
    config.batch_size = batch_size
    config.num_epoch = epoches
    config.use_mix = use_mix
    config.old_ratio = old_ratio
    config.starter_learning_rate = lr
    config.decay_steps = decay_steps
    config.decay_steps_inner = decay_steps_inner
    config.decay_rate = decay_rate
    config.restart = restart
    config.resnet = resnet
    return config


def train(
    cv_dim_list,
    neurons=[240, 120, 60, 30],
    numb_threads=8,
    resnet=True,
    use_mix=False,
    restart=False,
    batch_size=128,
    epoches=12000,
    lr=0.0008,
    decay_steps=120,
    decay_rate=0.96,
    old_ratio=7.0,
    decay_steps_inner=0,
    init_model=None
):
    config = set_conf(cv_dim_list,
                      neurons=neurons,
                      numb_threads=numb_threads,
                      resnet=resnet,
                      use_mix=use_mix,
                      restart=restart,
                      batch_size=batch_size,
                      epoches=epoches,
                      lr=lr,
                      decay_steps=decay_steps,
                      decay_rate=decay_rate,
                      old_ratio=old_ratio,
                      decay_steps_inner=decay_steps_inner)
    if init_model is not None:
        if config.restart:
            raise RuntimeError(
                "option --restart is conflicting with --init-model")
        config.graph_file = init_model
    reset_batch_size(config)
    print_conf(config, numb_threads)

    tf.reset_default_graph()
    tf_config = tf.ConfigProto(intra_op_parallelism_threads=numb_threads,
                               inter_op_parallelism_threads=2)
    with tf.Session(config=tf_config) as sess:
        print("Begin to optimize")
        reader = Reader(config)
        model = Model(config, sess)
        model.train(reader)


def get_parm():
    parser = argparse.ArgumentParser(
        description="*** Train a model. ***")
    parser.add_argument('-t', '--numb-threads', type=int, default=1,
                        help='the number of threads.')
    parser.add_argument('-c', '--cv-dim-list', type=int, default=[0, 0], nargs='+',
                        help='the list of cv dimension number')
    parser.add_argument('-r', '--restart', action='store_true',
                        help='if restart from a model.')
    parser.add_argument('-n', '--neurons', type=int, default=[240, 120, 60, 30], nargs='+',
                        help='the number of neurons in each hidden layer.')
    parser.add_argument('-b', '--batch-size', type=int, default=64,
                        help='the batch size.')
    parser.add_argument('-e', '--numb-epoches', type=int, default=3000,
                        help='the number of epoches.')
    parser.add_argument('-l', '--starter-lr', type=float, default=0.003,
                        help='the starter learning rate.')
    parser.add_argument('-m', '--use-mix', action='store_true',
                        help='mix the new data with old data.')
    parser.add_argument('-o', '--old-ratio', type=float, default=7,
                        help='the ratio of old data in a batch.')
    parser.add_argument('--decay-steps', type=int, default=10,
                        help='the decay steps.')
    parser.add_argument('--decay-steps-inner', type=int, default=0,
                        help='the inner loop decay steps. set 0 to disable inner lr loop')
    parser.add_argument('--decay-rate', type=float, default=0.96,
                        help='the decay rate.')
    parser.add_argument('--resnet', action='store_true',
                        help='try using resNet if two neighboring layers are of the same size.')
    parser.add_argument('-i', '--init-model', type=str,
                        help='use this graph to init the weights in the model.')
    args = parser.parse_args()
    return args


def main():
    args = get_parm()

    config = Config(args.cv_dim_list)
    config.n_neuron = args.neurons
    config.batch_size = args.batch_size
    config.num_epoch = args.numb_epoches
    config.use_mix = args.use_mix
    config.old_ratio = args.old_ratio
    config.starter_learning_rate = args.starter_lr
    config.decay_steps = args.decay_steps
    config.decay_steps_inner = args.decay_steps_inner
    config.decay_rate = args.decay_rate
    config.restart = args.restart
    config.resnet = args.resnet

    print("configuration:\n", config.__dict__)
    if args.init_model is not None:
        if config.restart:
            raise RuntimeError(
                "option --restart is conflicting with --init-model")
        config.graph_file = args.init_model
    reset_batch_size(config)
    print_conf(config, args.numb_threads)

    tf.reset_default_graph()
    tf_config = tf.ConfigProto(intra_op_parallelism_threads=args.numb_threads,
                               inter_op_parallelism_threads=2)
    with tf.Session(config=tf_config) as sess:
        print("Begin to optimize")
        reader = Reader(config)
        model = Model(config, sess)
        model.train(reader)


if __name__ == '__main__':
    main()
