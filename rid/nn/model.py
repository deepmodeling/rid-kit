import os
import time
import sys
import logging
try:
    import tensorflow.compat.v1 as tf
    tf.disable_v2_behavior()
except ImportError:
    import tensorflow as tf
import numpy as np
from tensorflow.python.ops import control_flow_ops
from tensorflow.python.training import moving_averages
from rid.constants import kbT, beta, N_grid, inverse_f_cvt

tf_precision = tf.float32


class Reader(object):
    def __init__(self, config):
        # copy from config
        self.data_path = config.data_path
        self.num_epoch = config.num_epoch
        self.use_mix = config.use_mix
        self.old_ratio = config.old_ratio
        assert self.old_ratio >= 0, "the old data ration sould be >= 0"
        self.batch_size = config.batch_size
        self.batch_size_old = int(
            self.batch_size * self.old_ratio / (1. + self.old_ratio))
        self.batch_size_new = self.batch_size - self.batch_size_old
        self.cv_dim = config.cv_dim
        self.drop_out_rate = config.drop_out_rate
        self.angular_mask = config.angular_mask

    def prepare(self):
        self.n_input = self.cv_dim
        self.index_count_all = 0
        self.index_count_new = 0
        self.index_count_old = 0
        if self.use_mix:
            tr_data_new = np.load(self.data_path)
            tr_data_new = np.reshape(tr_data_new, [-1, self.cv_dim * 2])
            tr_data_new[:, self.cv_dim:] *= inverse_f_cvt
            self.inputs_train_new = tr_data_new[:, :]
            tr_data_old = np.load(self.data_path)
            tr_data_old = np.reshape(tr_data_old, [-1, self.cv_dim * 2])
            tr_data_old[:, self.cv_dim:] *= inverse_f_cvt
            self.inputs_train_old = tr_data_old[:, :]
            self.train_size_new = self.inputs_train_new.shape[0]
            self.train_size_old = self.inputs_train_old.shape[0]
            if self.batch_size_new > self.train_size_new:
                self.batch_size_new = self.train_size_new
            if self.batch_size_old > self.train_size_old:
                self.batch_size_old = self.train_size_old
            self.batch_size = self.batch_size_old + self.batch_size_new
            print("# batch_size %d mixed by old %d new %d" %
                  (self.batch_size, self.batch_size_old, self.batch_size_new))
        else:
            tr_data_all = np.load(self.data_path)
            tr_data_all[:, self.cv_dim:] *= inverse_f_cvt
            self.inputs_train_all = tr_data_all[:, :]
            self.train_size_all = self.inputs_train_all.shape[0]
        # print(np.shape(self.inputs_train))

    def _sample_train_all(self):
        self.index_count_all += self.batch_size
        if self.index_count_all > self.train_size_all:
            # shuffle the data
            self.index_count_all = self.batch_size
            ind = np.random.choice(self.train_size_all,
                                   self.train_size_all, replace=False)
            self.inputs_train_all = self.inputs_train_all[ind, :]
        ind = np.arange(self.index_count_all -
                        self.batch_size, self.index_count_all)
        return self.inputs_train_all[ind, :]

    def _sample_train_mix(self, cat=True):
        self.index_count_new += self.batch_size_new
        if self.index_count_new > self.train_size_new:
            # shuffle the data
            self.index_count_new = self.batch_size_new
            ind = np.random.choice(self.train_size_new,
                                   self.train_size_new, replace=False)
            self.inputs_train_new = self.inputs_train_new[ind, :]

        self.index_count_old += self.batch_size_old
        if self.index_count_old > self.train_size_old:
            # shuffle the data
            self.index_count_old = self.batch_size_old
            ind = np.random.choice(self.train_size_old,
                                   self.train_size_old, replace=False)
            self.inputs_train_old = self.inputs_train_old[ind, :]

        ind_new = np.arange(self.index_count_new -
                            self.batch_size_new, self.index_count_new)
        ind_old = np.arange(self.index_count_old -
                            self.batch_size_old, self.index_count_old)
        if cat:
            return np.concatenate([self.inputs_train_new[ind_new, :], self.inputs_train_old[ind_old, :]], axis=0)
        else:
            return self.inputs_train_new[ind_new, :], self.inputs_train_old[ind_old, :]

    def sample_train(self, cat=True):
        if self.use_mix:
            return self._sample_train_mix(cat)
        else:
            return self._sample_train_all()

    def get_train_size(self):
        if self.use_mix:
            return int(self.train_size_new * (1. + self.old_ratio))
        else:
            return self.train_size_all

    def get_batch_size(self):
        return self.batch_size

    def get_data(self):
        if self.use_mix:
            return self.inputs_train_new, self.inputs_train_old
        else:
            return self.inputs_train_all


class Model(object):
    def __init__(self, config, sess):
        self.log_name = config.log_name
        self.sess = sess
        # copy from config
        self.data_path = config.data_path
        self.n_neuron = config.n_neuron
        self.useBN = config.useBN
        self.n_displayepoch = config.n_displayepoch
        self.starter_learning_rate = config.starter_learning_rate
        self.decay_steps = config.decay_steps
        self.decay_steps_inner = config.decay_steps_inner
        self.decay_rate = config.decay_rate
        self.display_in_training = config.display_in_training
        self.restart = config.restart
        self.resnet = config.resnet
        self.graph_file = config.graph_file
        self.cv_dim = int(config.cv_dim)
        self.angular_mask = np.array(config.angular_mask)
        self.angular_dim = np.sum(self.angular_mask, dtype=int)
        self.non_angular_dim = self.cv_dim - self.angular_dim
        self.angular_mask_boolean = (self.angular_mask == 1)
        self.non_angular_mask_boolean = (self.angular_mask == 0)
        if self.graph_file is not None:
            self.graph = self.load_graph(self.graph_file, 'load')
        else:
            self.graph = None

    def test_error(self, inputs_train):
        ret = self.sess.run([self.l2_loss, self.rel_error_k],
                            feed_dict={self.inputs_train: inputs_train,
                                       self.is_training: False,
                                       self.drop_out_rate: 0.0})
        error = np.sqrt(ret[0])
        error2 = np.mean(np.sqrt(ret[1]))
        return error, error2

    def test_error_mix(self, reader):
        data_new, data_old = reader.get_data()
        ret = self.sess.run([self.l2_loss, self.rel_error_k],
                            feed_dict={self.inputs_train: data_new,
                                       self.is_training: False,
                                       self.drop_out_rate: 0.0})
        error_new = np.sqrt(ret[0])
        error_new2 = np.mean(np.sqrt(ret[1]))
        ret = self.sess.run([self.l2_loss, self.rel_error_k],
                            feed_dict={self.inputs_train: data_old,
                                       self.is_training: False,
                                       self.drop_out_rate: 0.0})
        error_old = np.sqrt(ret[0])
        error_old2 = np.mean(np.sqrt(ret[1]))
        return error_new, error_new2, error_old, error_old2

    def train(self, reader):
        logger = logging.getLogger(__name__)
        logger.setLevel(level=os.environ.get("LOGLEVEL", "INFO").upper())
        handler = logging.FileHandler(self.log_name)
        formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s")
        handler.setFormatter(fmt=formatter)
        logger.addHandler(handler)
        
        reader.prepare()
        self.n_input = reader.n_input
        self.inputs_train = tf.placeholder(
            tf_precision, [None, self.n_input + self.cv_dim], name='inputs')
        self.is_training = tf.placeholder(tf.bool)
        self.drop_out_rate = tf.placeholder(tf.float32, name='drop_out_rate')  ### drop out ###
        self._extra_train_ops = []
        self.global_step = tf.get_variable('global_step', [],
                                           initializer=tf.constant_initializer(
                                               1),
                                           trainable=False, dtype=tf.int32)
        self.global_epoch = self.global_step * \
            reader.get_batch_size() // reader.get_train_size()
        if self.decay_steps_inner == 0:
            self.learning_rate = tf.train.exponential_decay(self.starter_learning_rate,
                                                            self.global_epoch,
                                                            self.decay_steps,
                                                            self.decay_rate,
                                                            staircase=True)
        else:
            self.global_epoch_inner = self.global_epoch % self.decay_steps
            self.lr_pref = tf.train.exponential_decay(1.0,
                                                      self.global_epoch,
                                                      self.decay_steps,
                                                      self.decay_rate,
                                                      staircase=False)
            self.learning_rate = tf.train.exponential_decay(self.starter_learning_rate,
                                                            self.global_epoch_inner,
                                                            self.decay_steps_inner,
                                                            self.decay_rate,
                                                            staircase=True)
            self.learning_rate *= self.lr_pref
        self.mv_decay = 1.0 - self.learning_rate/self.starter_learning_rate

        avg_input, scl_input = self.compute_statistic(reader)
        self.energy, self.l2_loss, self.rel_error_k\
            = self.build_force(self.inputs_train, suffix="test", reuse=False, shift=avg_input, scale=scl_input, graph=self.graph)

        # train operations
        trainable_variables = tf.trainable_variables()
        grads = tf.gradients(self.l2_loss, trainable_variables)
        optimizer = tf.train.AdamOptimizer(learning_rate=self.learning_rate)
        apply_op = optimizer.apply_gradients(zip(grads, trainable_variables),
                                             global_step=self.global_step, name='train_step')
        train_ops = [apply_op] + self._extra_train_ops
        self.train_op = tf.group(*train_ops)

        saver = tf.train.Saver()

        # initialization
        if (self.restart == False):
            sample_used = 0
            epoch_used = 0
            self.sess.run(tf.global_variables_initializer())
            print('# start training from scratch')
        elif self.graph is None:
            self.sess.run(tf.global_variables_initializer())
            saver.restore(self.sess, "old_model/model.ckpt")
            print("Model restored.")
            self.global_step = tf.assign(
                self.global_step, 0, name='global_step')
            cur_step = self.sess.run(self.global_step)
            sample_used = cur_step * reader.get_batch_size()
            epoch_used = sample_used // reader.get_train_size()

        start_time = time.time()

        inputs_train = reader.sample_train()
        if reader.use_mix:
            error, error2, error_old, error_old2 = self.test_error_mix(reader)
        else:
            error, error2 = self.test_error(inputs_train)
        current_lr = self.sess.run(tf.to_double(self.learning_rate))
        if self.display_in_training:
            if reader.use_mix:
                logger.info("epoch: %3u, ab_err_n: %.4e, rel_err_n: %.4e, ab_err_o: %.4e, rel_err_o: %.4e, lr: %.4e"
                      % (epoch_used, error, error2, error_old, error_old2, current_lr))
            else:
                logger.info("epoch: %3u, ab_err: %.4e, rel_err: %.4e, lr: %.4e" %
                      (epoch_used, error, error2, current_lr))

        while epoch_used < reader.num_epoch:
            # print('# doing training')
            inputs_train = reader.sample_train()
            self.sess.run([self.train_op],
                          feed_dict={self.inputs_train: inputs_train,
                                     self.is_training: True,
                                     self.drop_out_rate: reader.drop_out_rate})
            sample_used += reader.get_batch_size()
            # print(sample_used)
            if (sample_used // reader.get_train_size()) > epoch_used:
                epoch_used = sample_used // reader.get_train_size()
                if epoch_used % self.n_displayepoch == 0:
                    save_path = saver.save(
                        self.sess, os.getcwd() + "/" + "model.ckpt")
                    if reader.use_mix:
                        error, error2, error_old, error_old2 = self.test_error_mix(
                            reader)
                    else:
                        error, error2 = self.test_error(inputs_train)
                    current_lr = self.sess.run(
                        tf.to_double(self.learning_rate))
                    if self.display_in_training:
                        if reader.use_mix:
                            logger.info("epoch: %3u, ab_err_n: %.4e, rel_err_n: %.4e, ab_err_o: %.4e, rel_err_o: %.4e, lr: %.4e"
                                  % (epoch_used, error, error2, error_old, error_old2, current_lr))
                        else:
                            logger.info("epoch: %3u, ab_err: %.4e, rel_err: %.4e, lr: %.4e" % (
                                epoch_used, error, error2, current_lr))
                        sys.stdout.flush()
        end_time = time.time()
        logger.info("running time: %.3f s" % (end_time-start_time))

    def compute_statistic(self,
                          reader):
        max_scale = 3.
        if reader.use_mix:
            dnew, dold = reader.get_data()
        else:
            dold = reader.get_data()
        da = np.average(dold[:, 0:self.cv_dim], axis=0)
        ds = np.std(dold[:, 0:self.cv_dim], axis=0)
        # da[:self.cv_dih_dim] = 0.0
        da[self.angular_mask_boolean] = 0.0
        if all(ds != 0):
            ds = 1./(ds)
        # ds[:self.cv_dih_dim] = 1.0
        ds[self.angular_mask_boolean] = 1.0
        non_angular_cv = ds[self.non_angular_mask_boolean]
        non_angular_cv[non_angular_cv>max_scale] = max_scale
        ds[self.non_angular_mask_boolean] = non_angular_cv
        # for ii in range(self.cv_dih_dim, self.cv_dim):
        #     if ds[ii] > max_scale:
        #         ds[ii] = max_scale
        return da, ds

    def build_force(self,
                    inputs,
                    suffix,
                    shift=None,
                    scale=None,
                    reuse=None,
                    graph=None):
        cvs = tf.slice(inputs, [0, 0], [-1, self.cv_dim], name='cvs')
        if shift is not None:
            assert(scale is not None)
            t_shift = tf.get_variable('input_shift',
                                      [self.cv_dim],
                                      dtype=tf_precision,
                                      trainable=False,
                                      initializer=tf.constant_initializer(shift))
            t_scale = tf.get_variable('input_scale',
                                      [self.cv_dim],
                                      dtype=tf_precision,
                                      trainable=False,
                                      initializer=tf.constant_initializer(scale))
            # t_shift = tf.constant(shift, name='input_shift')
            # t_scale = tf.constant(scale, name='input_scale')
            cvs = (cvs - t_shift) * t_scale
        # angles = tf.slice(cvs, [0, 0], [-1, self.cv_dih_dim], name='angles')
        angles = tf.boolean_mask(
                cvs, self.angular_mask_boolean, axis=1, name='angles'
                )
        angles = tf.reshape(angles, [-1, self.angular_dim])
        dists = tf.boolean_mask(
                cvs, self.non_angular_mask_boolean, axis=1, name='dists'
                )
        dists = tf.reshape(dists, [-1, self.non_angular_dim])
        # dists = tf.slice(cvs, [0, self.cv_dih_dim],
        #                  [-1, self.cv_dist_dim], name='dists')
        forces_hat = tf.slice(
            inputs, [0, self.cv_dim], [-1, self.cv_dim], name='forces')
        inputs = tf.concat([tf.cos(angles), tf.sin(angles), dists], 1)
        if graph is not None:
            init_t = [graph.get_tensor_by_name('load/layer_0/matrix:0'),
                      graph.get_tensor_by_name('load/layer_0/bias:0')
                      ]
            with tf.Session(graph=self.graph) as g_sess:
                init = g_sess.run(init_t)
        else:
            init = None
        layer = self._one_layer(
            inputs, self.n_neuron[0], drop_out_rate=self.drop_out_rate, name='layer_0', reuse=reuse, init=init)
        for ii in range(1, len(self.n_neuron)):
            if graph is not None:
                init_t = [graph.get_tensor_by_name('load/layer_%s/matrix:0' % str(ii)),
                          graph.get_tensor_by_name(
                              'load/layer_%s/bias:0' % str(ii))
                          ]
                if self.resnet and self.n_neuron[ii] == self.n_neuron[ii-1]:
                    init_t += [graph.get_tensor_by_name(
                        'load/layer_%s/timestep:0' % str(ii))]
                with tf.Session(graph=self.graph) as g_sess:
                    init = g_sess.run(init_t)
            else:
                init = None
            if self.resnet and self.n_neuron[ii] == self.n_neuron[ii-1]:
                layer += self._one_layer(layer, self.n_neuron[ii], drop_out_rate=self.drop_out_rate, name='layer_'+str(
                    ii), reuse=reuse, with_timestep=True, init=init)
            else:
                layer = self._one_layer(
                    layer, self.n_neuron[ii], drop_out_rate=self.drop_out_rate, name='layer_'+str(ii), reuse=reuse, with_timestep=False, init=init)
        if graph is not None:
            init_t = graph.get_tensor_by_name('load/energy/matrix:0')
            with tf.Session(graph=self.graph) as g_sess:
                init = g_sess.run(init_t)
        else:
            init = None
        energy_ = self._final_layer(
            layer, 1, activation_fn=None, name='energy', reuse=reuse, init=init)
        energy = tf.identity(energy_, name='o_energy')
        energy_grad = tf.reshape(tf.stack(tf.gradients(energy, cvs)),
                              [-1, self.cv_dim], name='energy_grad')
        forces = tf.identity(-energy_grad, name='o_forces')
        force_dif = forces_hat - forces
        forces_norm = tf.reshape(tf.reduce_sum(
            forces * forces, axis=1), [-1, 1])
        forces_dif_norm = tf.reshape(tf.reduce_sum(
            force_dif * force_dif, axis=1), [-1, 1])
        l2_loss = tf.reduce_mean(forces_dif_norm, name='l2_loss')
        rel_error_k = forces_dif_norm / (1E-8 + forces_norm)
        return energy, l2_loss, rel_error_k

    def load_graph(self,
                   frozen_graph_filename,
                   prefix='load'):
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
                producer_op_list=None
            )
        # for op in graph.get_operations():
        #     print(op.name)
        return graph

    def _one_layer(self,
                   inputs,
                   outputs_size,
                   drop_out_rate,
                   activation_fn=tf.nn.tanh,
                   stddev=1.0,
                   bavg=0.0,
                   init=None,
                   name='linear',
                   reuse=None,
                   seed=None,
                   with_timestep=False):
        with tf.variable_scope(name, reuse=reuse):
            shape = inputs.get_shape().as_list()
            
            if init is not None:
                a_i_w = init[0]
                a_i_b = init[1]
                i_w_s = a_i_w.shape
                i_b_s = a_i_b.shape
                a_e_w = np.random.normal(
                    scale=stddev/np.sqrt(shape[1]+outputs_size), size=[shape[1], outputs_size])
                a_e_b = np.random.normal(
                    scale=stddev, loc=bavg, size=[outputs_size])
                a_e_w[:i_w_s[0], :i_w_s[1]] = a_i_w
                a_e_b[:i_b_s[0]] = a_i_b
                initer_w = tf.constant_initializer(a_e_w)
                initer_b = tf.constant_initializer(a_e_b)
            else:
                initer_w = tf.random_normal_initializer(
                    stddev=stddev/np.sqrt(shape[1]+outputs_size), seed=seed)
                initer_b = tf.random_normal_initializer(
                    stddev=stddev, mean=bavg, seed=seed)
            w = tf.get_variable('matrix',
                                [shape[1], outputs_size],
                                tf_precision,
                                initer_w)
            b = tf.get_variable('bias',
                                [outputs_size],
                                tf_precision,
                                initer_b)
            hidden = tf.matmul(inputs, w) + b
            if activation_fn != None and with_timestep:
                if init is not None:
                    a_i_t = init[2]
                    i_t_s = a_i_t.shape
                    a_e_t = np.random.normal(
                        scale=0.001, loc=0.1, size=[outputs_size])
                    a_e_t[0:i_t_s[0]] = a_i_t
                    initer_t = tf.constant_initializer(a_e_t)
                else:
                    initer_t = tf.random_normal_initializer(
                        stddev=0.001, mean=0.1, seed=seed)
                timestep = tf.get_variable('timestep',
                                           [outputs_size],
                                           tf_precision,
                                           initer_t)

        if activation_fn != None:
            if self.useBN:
                # None
                hidden_bn = self._batch_norm(
                    hidden, name=name+'_normalization', reuse=reuse)
                layer_out =  activation_fn(hidden_bn)
            else:
                if with_timestep:
                    layer_out = activation_fn(hidden) * timestep
                else:
                    layer_out = activation_fn(hidden)
        else:
            if self.useBN:
                # None
                layer_out = self._batch_norm(hidden, name=name+'_normalization', reuse=reuse)
            else:
                layer_out = hidden
        return tf.nn.dropout(layer_out, rate=drop_out_rate)


    def _final_layer(self,
                     inputs,
                     outputs_size,
                     activation_fn=tf.nn.tanh,
                     stddev=1.0,
                     bavg=0.0,
                     init=None,
                     name='linear',
                     reuse=None,
                     seed=None):
        with tf.variable_scope(name, reuse=reuse):
            shape = inputs.get_shape().as_list()
            if init is not None:
                a_i_w = init
                i_w_s = a_i_w.shape
                a_e_w = np.random.normal(
                    scale=stddev/np.sqrt(shape[1]+outputs_size), size=[shape[1], outputs_size])
                a_e_w[:i_w_s[0], :i_w_s[1]] = a_i_w
                initer_w = tf.constant_initializer(a_e_w)
            else:
                initer_w = tf.random_normal_initializer(
                    stddev=stddev/np.sqrt(shape[1]+outputs_size), seed=seed)
            w = tf.get_variable('matrix',
                                [shape[1], outputs_size],
                                tf_precision,
                                initer_w)
            hidden = tf.matmul(inputs, w)

        if activation_fn != None:
            if self.useBN:
                # None
                hidden_bn = self._batch_norm(
                    hidden, name=name+'_normalization', reuse=reuse)
                return activation_fn(hidden_bn)
            else:
                return activation_fn(hidden)
        else:
            if self.useBN:
                # None
                return self._batch_norm(hidden, name=name+'_normalization', reuse=reuse)
            else:
                return hidden

    def _batch_norm(self, x, name, reuse):
        """Batch normalization"""
        with tf.variable_scope(name, reuse=reuse):
            params_shape = [x.get_shape()[-1]]
            beta = tf.get_variable('beta', params_shape, tf_precision,
                                   initializer=tf.random_normal_initializer(0.0, stddev=0.1, dtype=tf_precision))
            gamma = tf.get_variable('gamma', params_shape, tf_precision,
                                    initializer=tf.random_uniform_initializer(0.1, 0.5, dtype=tf_precision))
        with tf.variable_scope(name+'moving', reuse=False):
            moving_mean = tf.get_variable('moving_mean', params_shape, tf_precision,
                                          initializer=tf.constant_initializer(
                                              0.0, tf_precision),
                                          trainable=False)
            moving_variance = tf.get_variable('moving_variance', params_shape, tf_precision,
                                              initializer=tf.constant_initializer(
                                                  1.0, tf_precision),
                                              trainable=False)
        # These ops will only be preformed when training
        mean, variance = tf.nn.moments(x, [0], name='moments')
        self._extra_train_ops.append(
            moving_averages.assign_moving_average(moving_mean, mean, self.mv_decay))
        self._extra_train_ops.append(moving_averages.assign_moving_average(
            moving_variance, variance, self.mv_decay))
        mean, variance = control_flow_ops.cond(self.is_training,
                                               lambda: (mean, variance),
                                               lambda: (moving_mean, moving_variance))
#       # elipson used to be 1e-5. Maybe 0.001 solves NaN problem in deeper net.
        y = tf.nn.batch_normalization(x, mean, variance, beta, gamma, 1e-6)
        y.set_shape(x.get_shape())
        return y
