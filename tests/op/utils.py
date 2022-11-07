import numpy as np

DATA_NEW = [[-1.579863999999999935e+00, 2.829998999999999931e+00, 6.681187999999999683e+00, -7.297216999999999842e+00]]
DATA_OLD = [[-1.467921000000000031e+00, 2.274875999999999898e+00, 6.123238999999999876e+00, 1.010345000000000049e+01], 
           [-2.258202999999999960e+00, -2.910984000000000016e+00, -3.138644000000000212e+00, -2.289113000000000042e+01]]
DATA_RAW = np.concatenate((DATA_OLD, DATA_NEW), axis=0)
DATA_EMPTY = []



def npy2txt(npy):
    data = np.load(npy)
    name = npy.split(".npy")[0]
    np.savetxt(name, data)


def txt2npy(txt):
    data = np.loadtxt(txt)
    name = str(txt)+".npy"
    np.save(name, data)