import os, sys
import logging
import argparse
try:
    import tensorflow.compat.v1 as tf
    tf.disable_v2_behavior()
except ImportError:
    import tensorflow as tf
import numpy as np
from rid.nn.model import Reader, Model


class Config(object):
    def __init__(self, cv_dim):
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
        self.data_path = './'
        self.restart = False
        self.resnet = False
        self.graph_file = None
        self.cv_dim = cv_dim
        self.display_in_training = True
        self.drop_out_rate = 0.5
        self.angular_mask = None


def reset_batch_size(config):
    tr_data = np.load(config.data_path)
    if tr_data.shape[0] < config.batch_size:
        config.batch_size = tr_data.shape[0]
        print("using new batch_size of %d" % config.batch_size)


def print_conf(config, nthreads):
    logger = logging.getLogger(__name__)
    logger.setLevel(level=os.environ.get("LOGLEVEL", "INFO").upper())
    handler = logging.FileHandler(config.log_name)
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s")
    handler.setFormatter(fmt=formatter)
    logger.addHandler(handler)
    
    logger.info("# restart           " + str(config.restart))
    logger.info("# num_threads       %d" % nthreads)
    logger.info("# neurons           " + str(config.n_neuron))
    logger.info("# batch size        " + str(config.batch_size))
    logger.info("# num_epoch         " + str(config.num_epoch))
    logger.info("# use_mix           " + str(config.use_mix))
    logger.info("# old_ratio         " + str(config.old_ratio))
    logger.info("# lr_0              " + str(config.starter_learning_rate))
    logger.info("# drop out rate     " + str(config.drop_out_rate))
    logger.info("# decay_steps       " + str(config.decay_steps))
    logger.info("# decay_steps_inner " + str(config.decay_steps_inner))
    logger.info("# decay_rate        " + str(config.decay_rate))
    logger.info("# resnet            " + str(config.resnet))
    logger.info("# graph_file        " + str(config.graph_file))


def set_conf(cv_dim,
             angular_mask,
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
             drop_out_rate=0.5,
             data_path="./",
             log_name = "log"):
    config = Config(cv_dim)
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
    config.drop_out_rate = drop_out_rate
    config.angular_mask = angular_mask
    config.data_path = data_path
    config.log_name = log_name
    return config


def train(
        cv_dim,
        angular_mask,
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
        init_model=None,
        drop_out_rate=0.5,
        data_path="./",
        log_name = "log"
    ):
    config = set_conf(cv_dim,
                      angular_mask=angular_mask,
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
                      decay_steps_inner=decay_steps_inner,
                      drop_out_rate = drop_out_rate,
                      data_path = data_path,
                      log_name = log_name)
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
    parser.add_argument('-d', '--drop-out-rate', type=float, default=0.5,
                        help='drop out rate')
    args = parser.parse_args()
    return args

