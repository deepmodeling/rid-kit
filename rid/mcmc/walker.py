import numpy as np
import tensorflow as tf

# kinetic enery in eV
kbT = (8.617343E-5) * 300 
beta = 1.0 / kbT

# ev to kj/mol
f_cvt = 96.485

def distance(cv1,cv2):
    d = np.sum((cv1 - cv2)**2,axis=1)
    return d

def s(cv, lm, xlist):
    s = 0
    n = 1e-10
    for i in range(1,10):
       s +=  i*np.exp(-lm*distance(cv,xlist[i-1]))
       n += np.exp(-lm*distance(cv,xlist[i-1]))
    return s/n

def z(cv, lm, xlist):
    z = 0
    for i in range(9):
        z += np.exp(-lm*distance(cv,xlist[i]))
    p = np.ones(shape=z.shape)*1e-10
    z = np.log(np.max((z,p),axis=0))
    return -z/lm

class Walker(object):
    def __init__(self, fd, nw, sess, type, cv_lower, cv_upper):
        self.type = type
        self._full_dim = fd
        self._num_walker = nw
        self._move_scale = 0.5
        self._sess = sess
        self._shape = (self._num_walker, self._full_dim)
        self._cv_range = cv_upper
        self._cv_shift = cv_lower
        # absolute coordinate
        self._position = np.random.uniform(size=self._shape)*(np.array(self._cv_range) - np.array(self._cv_shift)) + np.array(self._cv_shift)
        self._energy = np.zeros([self._num_walker])
        self._force = np.zeros([self._num_walker, self._full_dim])
        self._sample_step = 20
        self._acp_ratio_lb = 0.15
        self._acp_ratio_ub = 0.75
        self.max_scale = 5
        self.min_scale = 0.01
        self.inc_scale_fac = 1.25

    def sample(self, compute_ef, inter_step=1):
        acp_ratio = []
        self._energy, self._force = compute_ef(self._sess, self._position)
        for _ in range(inter_step):
            if self.type == "dih":
                position_new = np.mod(self._position + np.random.normal(scale=self._move_scale,size=self._shape), 2*np.pi)
            elif self.type == "dis":
                position_new = self._position + np.random.normal(scale=self._move_scale*0.08, size=self._shape)
                indices = np.where(position_new < 0.1)
                position_new[indices]=0.13 + np.random.normal(scale=0.01, size=indices[0].shape)
                indices = np.where(position_new > 9.9)
                position_new[indices]=9.87 + np.random.normal(scale=0.01, size=indices[0].shape)
            else:
                raise ValueError("Undefined cv type, only support 'dih' and 'dis' type")
            energy_new, force_new = compute_ef(self._sess, position_new)
            # in case of overflow
            prob_ratio = np.exp(np.minimum(- beta * (energy_new - self._energy), 0))
            idx = np.random.uniform(size=self._num_walker) < np.reshape(prob_ratio, [-1])

            self._position[idx, :] = position_new[idx, :]
            self._energy[idx] = energy_new[idx]
            self._force[idx] = force_new[idx]
            acp_ratio.append(np.mean(idx))
            # elif self.type == "dis":
            #     for ii in range(self._num_walker):
            #         if (position_new[ii,:] > self._cv_shift).all() and (position_new[ii,:] < self._cv_range).all() and idx[ii]:
            #             self._position[ii, :] = position_new[ii, :]
            #             self._energy[ii] = energy_new[ii]
            #             self._force[ii] = force_new[ii]
            #             acp_ratio.append(np.mean(idx))
            # else:
            #     raise ValueError("Undefined cv type, only support 'dih' and 'dis' type")
        acp_ratio = np.mean(acp_ratio)
        if acp_ratio > self._acp_ratio_ub:
            # move_scale is too small
            self._move_scale = np.min((self._move_scale*self.inc_scale_fac,self.max_scale),axis=0)
            print(
                "Increase move_scale to %s due to high acceptance ratio: %f" % (
                    self._move_scale, acp_ratio))
            # print(self._position[:5, :, :])
        elif acp_ratio < self._acp_ratio_lb:
            # move_scale is too large
            self._move_scale = np.max((self._move_scale/self.inc_scale_fac,self.min_scale),axis=0)
            print(
                "Decrease move_scale to %s due to low acceptance ratio: %f" % (
                    self._move_scale, acp_ratio))
        return self._position, self._energy, self._force

# project on all 1d points
def my_hist1d(pp, xx, delta, fd):
    my_hist = np.zeros((fd, len(xx)))
    for ii in range(pp.shape[0]):   ###trj_num
        for jj in range(fd):        ###cv_num
            my_hist[jj, np.int(pp[ii,jj]//delta)] += 1
    my_hist /= (pp.shape[0] * delta)
    return my_hist

# project on specific 2d point
def my_hist2d(pp, xx, yy, delta, cv1, cv2):
    my_hist = np.zeros((1, len(xx), len(yy)))
    for ii in range(pp.shape[0]):
        my_hist[0, np.int(pp[ii,cv1]//delta), np.int(pp[ii,cv2]//delta)] += 1
    my_hist /= (pp.shape[0] * delta * delta)
    return my_hist

# project on all path CV 1d
def my_hist1d_path(pp, xx, delta, lm, xlist, proj_index):
    my_hist = np.zeros((2, len(xx)))
    cv = pp[:,proj_index]
    slist = s(cv, lm, xlist)
    zlist = z(cv, lm, xlist)
    for ii in range(cv.shape[0]):
        if 0 < slist[ii] < 10 and -10 < zlist[ii] < 0:
            my_hist[0, np.int((slist[ii])//delta)] += 1
            my_hist[1, np.int((zlist[ii]+10)//delta)] += 1
    my_hist /= (pp.shape[0] * delta)
    return my_hist

# project on path CV 2d
def my_hist2d_path(pp, xx, yy, delta, lm, xlist, proj_index):
    my_hist = np.zeros((1, len(xx), len(yy)))
    cv = pp[:,proj_index]
    slist = s(cv, lm, xlist)
    zlist = z(cv, lm, xlist)
    for ii in range(cv.shape[0]):
        if 0 < slist[ii] < 10 and -10 < zlist[ii] < 0:
            my_hist[0, np.int((slist[ii])//delta), np.int((zlist[ii]+10)//delta)] += 1
    my_hist /= (pp.shape[0] * delta * delta)
    return my_hist