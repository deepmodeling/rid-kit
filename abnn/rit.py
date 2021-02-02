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
import time
# global consts
from lib.modeling import cv_dim
from lib.modeling import enhc_name
from lib.modeling import enhc_out_conf
from lib.modeling import enhc_out_angle
from lib.modeling import train_name
from lib.modeling import mol_name
from lib.modeling import mol_files
# utils
from lib.modeling import make_iter_name
from lib.modeling import make_walker_name
from lib.modeling import record_iter
from lib.modeling import log_iter
from lib.modeling import log_task
from lib.modeling import replace
from lib.modeling import make_grompp_enhc
from lib.modeling import copy_file_list
from lib.modeling import create_path
from lib.modeling import cmd_append_log
# tasks
from lib.modeling import make_res
from lib.modeling import run_res
from lib.modeling import post_res
from lib.modeling import make_train
from lib.modeling import run_train
# machine
import lib.MachineLocal as MachineLocal
import lib.MachineSlurm as MachineSlurm
from lib.machine_exec import exec_hosts
from lib.machine_exec import exec_hosts_batch

from template.tools.model_dev import compute_model_dev

exec_machine = MachineLocal
max_tasks = 1000000

temp_name=enhc_name
temp_out_conf=enhc_out_conf
temp_out_angle=enhc_out_angle
temp_files=["plumed.dat", "test.std.py"]
temp_plm="plumed.dat"
temp_out_plm="plm.out"
            
def make_temp (iter_index, 
               json_file, 
               graph_files) :
    graph_files.sort()
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    numb_walkers = jdata["numb_walkers"]
    template_dir = jdata["template_dir"]    
    nsteps = jdata["temp_nsteps"]
    frame_freq = jdata["temp_frame_freq"]
    start_temp = jdata["start_temp"]

    iter_name = make_iter_name (iter_index)
    work_path = iter_name + "/" + temp_name + "/"  
    mol_path = template_dir + "/" + mol_name + "/"
    temp_path = template_dir + "/" + temp_name + "/"
    conf_list = glob.glob(mol_path + "conf*gro")
    conf_list.sort()
    assert (len(conf_list) >= numb_walkers), "not enough conf files in mol dir %s" % mol_path

    create_path (work_path)

    for walker_idx in range (numb_walkers) :
        walker_path = work_path + make_walker_name(walker_idx) + "/"
        create_path (walker_path)
        # copy md ifles
        for ii in mol_files :
            if os.path.exists (walker_path + ii) :
                os.remove (walker_path + ii)
            shutil.copy (mol_path + ii, walker_path)
        # copy conf file
        conf_file = conf_list[walker_idx]
        if os.path.exists (walker_path + "conf.gro") :
            os.remove (walker_path + "conf.gro")
        shutil.copy (conf_file, walker_path + "conf.gro")
        # if have prev confout.gro, use as init conf
        if (iter_index > 0) :
            prev_temp_path = make_iter_name (iter_index-1) + "/" + temp_name + "/" + make_walker_name(walker_idx) + "/"
            prev_temp_path = os.path.abspath (prev_temp_path) + "/"
            if os.path.isfile (prev_temp_path + "confout.gro") :
                os.remove (walker_path + "conf.gro")
                os.symlink (prev_temp_path + "confout.gro", walker_path + "conf.gro")
            log_task ("use conf of iter " + make_iter_name(iter_index-1) + " walker " + make_walker_name(walker_idx) )
        # copy temp file
        for ii in temp_files :
            if os.path.exists (walker_path + ii) :
                os.remove (walker_path + ii)
            shutil.copy (temp_path + ii, walker_path)
        # copy graph files
        for ii in graph_files :
            file_name = os.path.basename(ii)
            abs_path = os.path.abspath(ii)        
            if os.path.exists (walker_path + file_name) :
                os.remove (walker_path + file_name)
            os.symlink (abs_path, walker_path + file_name)            
        # config MD
        mol_conf_file = walker_path + "grompp.mdp"
        make_grompp_enhc (mol_conf_file, nsteps, frame_freq)    
        # config plumed
        if iter_index == 0: 
            cur_temp = start_temp
        else :            
            cur_temp = np.loadtxt (os.path.join(prev_temp_path, "next.temp"))
        log_task (("use temp of %f") % (cur_temp))
        log_task (("length of traj %d") % (nsteps))
        np.savetxt (os.path.join(walker_path, 'cur.temp'), [cur_temp])
        plm_conf = walker_path + temp_plm
        replace(plm_conf, "TEMP=[^ ]* ", ("TEMP=%s " % cur_temp))
        replace(plm_conf, "STRIDE=[^ ]* ", ("STRIDE=%d " % frame_freq))
        replace(plm_conf, "FILE=[^ ]* ", ("FILE=%s " % temp_out_plm))

def run_temp (iter_index,
              json_file) :
    iter_name = make_iter_name (iter_index)
    work_path = iter_name + "/" + temp_name + "/"  

    fp = open (json_file, 'r')
    jdata = json.load (fp)
    gmx_prep = jdata["gmx_prep"]
    gmx_run = jdata["gmx_run"]
    temp_thread = jdata["temp_thread"]
    gmx_run = gmx_run + (" -nt %d " % temp_thread)
    gmx_prep_log = "gmx_grompp.log"
    gmx_run_log = "gmx_mdrun.log"
    gmx_run = gmx_run + " -plumed " + temp_plm
    gmx_prep_cmd = cmd_append_log (gmx_prep, gmx_prep_log)
    gmx_run_cmd = cmd_append_log (gmx_run, gmx_run_log)
    numb_walkers = jdata["numb_walkers"]

    all_task = glob.glob(work_path + "/[0-9]*[0-9]")
    all_task.sort()

    global exec_machine
    exec_hosts(MachineLocal, gmx_prep_cmd, 1, all_task, None)
    if len(all_task) == 1 :
        exec_hosts(MachineLocal, gmx_run_cmd, temp_thread, all_task, None)
    else :
        exec_hosts_batch(exec_machine, gmx_run_cmd, temp_thread, all_task, None)

def make_temp_decision (model_dev,
                        dev_trust_lvl,
                        temp_incr_lvl,
                        temp_decr_lvl) :
    nframes = model_dev.shape[0]
    untrust_ratio = np.sum(model_dev > dev_trust_lvl) / np.double(nframes)
    if untrust_ratio > temp_decr_lvl :
        return 1        # too many untrusted, decrease temp
    elif untrust_ratio < temp_incr_lvl :
        return 2        # too few untrusted, increase temp
    else :
        return 0        # keep temp

def make_next_temp (graph_list,
                    cur_temp, 
                    xx,
                    temp_test_intvl, 
                    dev_trust_lvl,
                    temp_incr_ratio, 
                    temp_decr_ratio,
                    temp_incr_lvl,
                    temp_decr_lvl, 
                    temp_min,
                    temp_max) :
    nframes = xx.shape[0]
    model_dev = compute_model_dev (graph_list, xx)
    assert(model_dev.shape[0] == nframes)

    decision_list = []
    for ii in range(0, nframes-1, temp_test_intvl) :
        sample_dev = model_dev[ii:ii+temp_test_intvl]        
        d = make_temp_decision(sample_dev, dev_trust_lvl, temp_incr_lvl, temp_decr_lvl)
        decision_list.append(d)
    decision_list = np.array(decision_list)
    
    if any(decision_list == 1) :
        next_temp = cur_temp * temp_decr_ratio
    elif all(decision_list == 2) :
        next_temp = cur_temp * temp_incr_ratio
    else :
        next_temp = cur_temp

    # cap temperature
    if next_temp < temp_min :
        next_temp = temp_min
    elif next_temp > temp_max :
        next_temp = temp_max
    
    return next_temp, decision_list, model_dev

def post_temp (iter_index, 
               json_file) :
    iter_name = make_iter_name (iter_index)
    work_path = iter_name + "/" + temp_name + "/"  
    
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    gmx_split = jdata["gmx_split_traj"]
    gmx_split_log = "gmx_split.log"
    gmx_split_cmd = cmd_append_log (gmx_split, gmx_split_log)
    temp_test_intvl = jdata['temp_test_intvl']
    dev_trust_lvl = jdata['dev_trust_lvl']
    temp_incr_ratio = jdata['temp_incr_ratio']
    temp_decr_ratio = jdata['temp_decr_ratio']
    temp_incr_lvl = jdata['temp_incr_lvl']
    temp_decr_lvl = jdata['temp_decr_lvl']
    phys_temp = jdata['phys_temp']
    max_temp = jdata['max_temp']
    
    all_task = glob.glob(work_path + "/[0-9]*[0-9]")
    all_task.sort()

    cwd = os.getcwd()
    numb_walkers = jdata["numb_walkers"]
    for ii in range(numb_walkers) :
        walker_path = work_path + make_walker_name(ii) + "/"
        os.chdir(walker_path)        
        if os.path.isdir ("confs") : 
            shutil.rmtree ("confs")
        os.makedirs ("confs")
        os.chdir(cwd)

    global exec_machine
    exec_hosts (MachineLocal, gmx_split_cmd, 1, all_task, None)

    walker_stop_flags = []
    for ii in range(numb_walkers) :
        walker_path = work_path + make_walker_name(ii) + "/"
        angles = np.loadtxt (walker_path + temp_out_plm)
        angles = angles[:,1:]
        np.savetxt (walker_path + temp_out_angle, angles, fmt = "%.6f")
        graph_files = glob.glob (walker_path + "/*.pb")
        cur_temp = np.loadtxt(os.path.join(walker_path, 'cur.temp'))
        if len (graph_files) != 0 :
            next_temp, dec_list, model_dev\
                = make_next_temp(graph_files, 
                                 cur_temp, 
                                 angles, 
                                 temp_test_intvl,
                                 dev_trust_lvl,
                                 temp_incr_ratio,
                                 temp_decr_ratio,
                                 temp_incr_lvl,
                                 temp_decr_lvl, 
                                 phys_temp, 
                                 max_temp)
            model_dev = np.concatenate((angles, np.reshape(model_dev, [-1,1])), axis = 1)
            np.savetxt(os.path.join(walker_path,'model.devi'), model_dev)
            np.savetxt(os.path.join(walker_path,'dec_list'), dec_list)
        else :        
            next_temp = cur_temp
        np.savetxt(os.path.join(walker_path,'next.temp'), [next_temp])

        if cur_temp == max_temp and next_temp == max_temp :
            walker_stop_flags.append(True)
        else :
            walker_stop_flags.append(False)

    # return if the iteration continues
    if all(walker_stop_flags) :
        return False
    else :
        return True

def run_iter (json_file, init_model) :
    prev_model = init_model
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    numb_iter = jdata["numb_iter"]
    numb_task = 8
    record = "record.rit"

    iter_rec = [0, -1]
    if os.path.isfile (record) :
        with open (record) as frec :
            for line in frec : 
                iter_rec = [int(x) for x in line.split()]
        logging.info ("continue from iter %03d task %02d" % (iter_rec[0], iter_rec[1]))

    global exec_machine

    for ii in range (numb_iter) :
        if ii > 0 :
            prev_model = glob.glob (make_iter_name(ii-1) + "/" + train_name + "/*pb")
        for jj in range (numb_task) :
            if ii * max_tasks + jj <= iter_rec[0] * max_tasks + iter_rec[1] : 
                continue
            if   jj == 0 :
                log_iter ("make_temp", ii, jj)
                make_temp (ii, json_file, prev_model)            
            elif jj == 1 :
                log_iter ("run_temp", ii, jj)
                run_temp  (ii, json_file)
            elif jj == 2 :
                log_iter ("post_temp", ii, jj)
                cont = post_temp (ii, json_file)
                if not cont : 
                    log_iter ("no more conf needed", ii, jj)
                    return
            elif jj == 3 :
                log_iter ("make_res", ii, jj)
                cont = make_res (ii, json_file)
            elif jj == 4 :
                log_iter ("run_res", ii, jj)
                run_res   (ii, json_file, exec_machine)
            elif jj == 5 :
                log_iter ("post_res", ii, jj)
                post_res  (ii, json_file)
            elif jj == 6 :
                log_iter ("make_train", ii, jj)
                make_train(ii, json_file)
            elif jj == 7 :
                log_iter ("run_train", ii, jj)
                run_train (ii, json_file, exec_machine)
            else :
                raise RuntimeError ("unknow task %d, something wrong" % jj)

            record_iter (record, ii, jj)


def _main () :
    parser = argparse.ArgumentParser()
    parser.add_argument("JSON", type=str, 
                        help="The json parameter")
    parser.add_argument("-m", "--models", default=[], nargs = '*', type=str, 
                        help="The init guessed model")    
    parser.add_argument("--machine", type=str, 
                        help="The machine settings")        
    args = parser.parse_args()

    logging.basicConfig (level=logging.INFO, format='%(asctime)s %(message)s')
    # logging.basicConfig (filename="compute_string.log", filemode="a", level=logging.INFO, format='%(asctime)s %(message)s')

    machine_type = "local"    
    gmxrc = None
    vcores = None
    if args.machine != None :
        fp = open (args.machine, 'r')
        jdata = json.load (fp)
        machine_type = jdata["machine_type"]
        gmxrc = jdata["gmxrc"]
        vcores = jdata["virtual_cores"]

    global exec_machine
    if   machine_type == "local" :
        exec_machine = MachineLocal
    elif machine_type == "slurm" :
        exec_machine = MachineSlurm

    if vcores != None:
        exec_machine.has_virtual_cores(vcores)
    if gmxrc != None:
        exec_machine.add_source_file(gmxrc)

    logging.info ("start running")
    run_iter (args.JSON, args.models)
    logging.info ("finished!")

if __name__ == '__main__':
    _main()
