#!/usr/bin/env python3

import os
import re
import shutil
import json
import argparse
import numpy as np
import subprocess as sp
import glob
import logging
import librun as libr
import libstring as libs

import MachineLocal, MachineSlurm

from strm import *

def run_iter (json_file) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    numb_iter = jdata["numb_iter"]
    tol = jdata["tol_diff"]
    init_method = jdata["init_method"]
    record = "record.strm"
    numb_task = 5
    max_tasks = 1000000

    iter_rec = [0, -1]
    if os.path.isfile (record) :
        with open (record) as frec :
            for line in frec : 
                iter_rec = [int(x) for x in line.split()]
        logging.info ("continue from iter %03d task %02d" % (iter_rec[0], iter_rec[1]))

    if iter_rec[0] * max_tasks + iter_rec[1] >= 0 :
        return
    ii = 0
    jj = 0
    if init_method == "confs_linear" :
        libr.log_iter ("init_string with confs_linear", ii, jj)
        init_string_confs_linear(json_file)
        libr.log_iter ("equi_string with confs_linear", ii, jj)
        equi_0_serial (json_file)
    elif init_method == "resample_old" :
        libr.log_iter ("init_string with resample_old", ii, jj)
        init_string_resample_old(json_file)
        libr.log_iter ("equi_string with resample_old", ii, jj)
        equi_old_string_parallel (json_file)
    else :
        raise RuntimeError("unknown init method: " + init_method)
    libr.record_iter (record, ii, jj)            

def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument("JSON", type=str, 
                        help="The json parameter")
    args = parser.parse_args()

    logging.basicConfig (level=logging.INFO, format='%(asctime)s %(message)s')

    logging.info ("start running")
    run_iter (args.JSON)
    logging.info ("finished!")    

if __name__ == '__main__':
    _main()
