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
import lib.MachineLocal as MachineLocal
import lib.MachineSlurm as MachineSlurm
# machine
import lib.MachineLocal as MachineLocal
import lib.MachineSlurm as MachineSlurm
from lib.machine_exec import exec_hosts
from lib.machine_exec import exec_hosts_batch
from lib.batch_exec import exec_batch
from lib.batch_exec import exec_batch_group

from template.tools.cluster_cv import sel_from_cluster

cv_dim = 4
cv_dih_dim = 0

shell_clustering = False

iter_format = "%06d"
walker_format = "%03d"
task_format = "%02d"
log_iter_head = "iter " + iter_format + " task " + task_format + ": "

enhc_name="00.enhcMD"
enhc_out_conf="confs/"
enhc_out_angle="angle.rad.out"

mol_name="mol"
mol_files=["grompp.mdp", "topol.top"]

res_name="01.resMD"
res_files=["cmpf.sh", "cmpf.py", "general_mkres.sh", "plumed.res.templ", "tools"]
res_plm="plumed.res.dat"

train_name="02.train"
train_files=["model.py", "main.py", "freeze.py"]

def make_iter_name (iter_index) :
    return "iter." + (iter_format % iter_index)

def make_walker_name (walker_index) :
    return (walker_format % walker_index)

def record_iter (record, ii, jj) :
    with open (record, "a") as frec :
        frec.write ("%d %d\n" % (ii, jj))        

def log_iter (task, ii, jj) :
    logging.info ((log_iter_head + "%s") % (ii, jj, task))

def repeat_to_length(string_to_expand, length):
    ret = ""
    for ii in range (length) : 
        ret += string_to_expand
    return ret

def log_task (message) :
    header = repeat_to_length (" ", len(log_iter_head % (0, 0)))
    logging.info (header + message)

def clean_files(files) :
    for ii in files :
        jlist = glob.glob(ii)
        for jj in jlist:
            if os.path.isdir(jj) :
                shutil.rmtree(jj)
            else :
                os.remove(jj)            
    
def replace (file_name, pattern, subst) :
    file_handel = open (file_name, 'r')
    file_string = file_handel.read ()
    file_handel.close ()
    file_string = ( re.sub (pattern, subst, file_string) )
    file_handel = open (file_name, 'w')
    file_handel.write (file_string)
    file_handel.close ()

def make_grompp_enhc (gro_file, nsteps, frame_freq) :
    replace (gro_file, "nsteps.*=.*", "nsteps = %d" % nsteps)
    replace (gro_file, "nstxout.*=.*", "nstxout = %d" % frame_freq)
    replace (gro_file, "nstvout.*=.*", "nstvout = %d" % frame_freq)
    replace (gro_file, "nstfout.*=.*", "nstfout = %d" % frame_freq)
    replace (gro_file, "nstxtcout.*=.*", "nstxtcout = %d" % frame_freq)
    replace (gro_file, "nstenergy.*=.*", "nstenergy = %d" % frame_freq)    

def make_grompp_res (gro_file, nsteps, frame_freq) :
    replace (gro_file, "nsteps.*=.*", "nsteps = %d" % nsteps)
    replace (gro_file, "nstxout.*=.*", "nstxout = %d" % 0)
    replace (gro_file, "nstvout.*=.*", "nstvout = %d" % 0)
    replace (gro_file, "nstfout.*=.*", "nstfout = %d" % 0)
    replace (gro_file, "nstenergy.*=.*", "nstenergy = %d" % 0)

def copy_file_list (file_list, from_path, to_path) :
    for jj in file_list : 
        if os.path.isfile(from_path + jj) :
            shutil.copy (from_path + jj, to_path)
        elif os.path.isdir(from_path + jj) :
            cwd=os.getcwd()
            os.chdir(from_path+jj)
            files = glob.glob("*")
            os.chdir(cwd)
            os.makedirs(to_path+jj)
            for ff in files :
                shutil.copy(from_path+jj+'/'+ff, to_path+jj+'/'+ff)
            #    print(from_path+jj+'/'+ff, to_path+jj+'/'+ff)
            # if os.path.exists(to_path+jj) :
            #     print('error!!!!!!!!!!!!!!!!')
            # shutil.copytree (from_path + jj, to_path + jj)

def create_path (path) :
    path += '/'
    if os.path.isdir(path) : 
        dirname = os.path.dirname(path)
        counter = 0
        while True :
            bk_dirname = dirname + ".bk%03d" % counter
            if not os.path.isdir(bk_dirname) : 
                shutil.move (dirname, bk_dirname) 
                break
            counter += 1
    os.makedirs (path)

def cmd_append_log (cmd,
                    log_file) :
    ret = cmd
    ret = ret + " 1> " + log_file
    ret = ret + " 2> " + log_file
    return ret

def make_res (iter_index, 
              json_file) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    numb_walkers = jdata["numb_walkers"]
    template_dir = jdata["template_dir"]    
    bias_nsteps = jdata["bias_nsteps"]
    bias_frame_freq = jdata["bias_frame_freq"]
    nsteps = jdata["res_nsteps"]
    frame_freq = jdata["res_frame_freq"]
    sel_threshold = jdata["sel_threshold"]
    max_sel = jdata["max_sel"]
    cluster_threshold = jdata["cluster_threshold"]
    init_numb_cluster=[16,26]
    base_path = os.getcwd() + "/"
    iter_name = make_iter_name (iter_index)
    enhc_path = iter_name + "/" + enhc_name + "/" 
    os.chdir (enhc_path)
    enhc_path = os.getcwd() + "/"
    os.chdir (base_path)
    templ_mol_path = template_dir + "/" + mol_name + "/"
    templ_res_path = template_dir + "/" + res_name + "/"
    res_path = iter_name + "/" + res_name + "/"
    create_path (res_path) 

    ret_list = [True for ii in range(numb_walkers)]
    
    # sel angles
    ## check if we have graph in enhc
    for walker_idx in range(numb_walkers) :
        walker_path = enhc_path + walker_format % walker_idx + "/"
        graph_files = glob.glob (walker_path + "/*.pb")
        if len (graph_files) != 0 :
            cluster_threshold = np.loadtxt("cluster_threshold.dat")
            os.chdir (walker_path)
            sel_cmd = "python3 test.std.py -m *.pb -t %f -d %s --output sel.out --output-angle sel.angle.out" % (sel_threshold, enhc_out_angle)
            sel_cmd = cmd_append_log (sel_cmd, "sel.log")
            log_task ("select with threshold %f" % sel_threshold)
            log_task (sel_cmd)
            sp.check_call (sel_cmd, shell = True)
            os.chdir (base_path)
            sel_idx = []
            sel_angles = np.array([])
            with open (walker_path + "sel.out") as fp :
                for line in fp : 
                    sel_idx += [int(x) for x in line.split()]
            if len(sel_idx) != 0 :
                sel_angles = np.reshape(np.loadtxt(walker_path + 'sel.angle.out'), [-1,cv_dim])
            elif len(sel_idx) == 0 :
                np.savetxt (walker_path + 'num_of_cluster.dat', [0], fmt = '%d')
                np.savetxt (walker_path + 'cls.sel.out', [], fmt = '%d')
                continue

        else :
            cluster_threshold = jdata["cluster_threshold"]
            sel_idx = range (len(glob.glob (walker_path + enhc_out_conf + "conf*gro")))
            sel_angles = np.loadtxt (walker_path + enhc_out_angle)
            sel_angles = np.reshape (sel_angles, [-1, cv_dim])            
            np.savetxt(walker_path + 'sel.out', sel_idx, fmt = '%d')
            np.savetxt(walker_path + 'sel.angle.out', sel_angles, fmt = '%.6f')
            if walker_idx ==0:
                cls_sel = sel_from_cluster (sel_angles, cluster_threshold)
                test_numb_cluster=len(set(cls_sel))
                print(test_numb_cluster)
                for test_iter in range(500):
                    if test_numb_cluster < init_numb_cluster[0]:
                        cluster_threshold=cluster_threshold*0.95
                        cls_sel = sel_from_cluster (sel_angles, cluster_threshold)
                        test_numb_cluster=len(set(cls_sel))
                    elif test_numb_cluster > init_numb_cluster[1]:
                        cluster_threshold=cluster_threshold*1.05
                        cls_sel = sel_from_cluster (sel_angles, cluster_threshold)
                        test_numb_cluster=len(set(cls_sel))
                    else:
                        np.savetxt (walker_path + 'cluster_threshold.dat', [cluster_threshold], fmt = '%f')
                        np.savetxt ('cluster_threshold.dat', [cluster_threshold], fmt = '%f')
                        break
            else:
                cluster_threshold = np.loadtxt("cluster_threshold.dat")
                
        conf_start = 0
        conf_every = 1

        sel_idx = np.array (sel_idx, dtype = np.int)
        assert (len(sel_idx) == sel_angles.shape[0])
        if shell_clustering and len(sel_idx)>1:
            cmd_sel_from_cluster = (base_path +
                                    "template/tools/cluster_cv.py -i %s -c %s -t %f --output-idx %s  --output-cv %s" 
                                    % ( walker_path + 'sel.out',
                                        walker_path + 'sel.angle.out',
                                        cluster_threshold,
                                        walker_path + 'cls.sel.out',
                                        walker_path + 'cls.sel.angle.out'
                                    )
            )
            sp.check_call(cmd_sel_from_cluster, shell = True)
            sel_idx = np.loadtxt(walker_path + 'cls.sel.out', dtype = np.int)
        elif shell_clustering == False and len(sel_idx)>1 :
            cls_sel = sel_from_cluster (sel_angles, cluster_threshold)
##############################################################################################
            np.savetxt (walker_path + 'num_of_cluster.dat', [len(set(cls_sel))], fmt = '%d')
            np.savetxt (walker_path + 'cluster_threshold.dat', [cluster_threshold], fmt = '%f')
            if len(cls_sel)>max_sel:
                cls_sel=cls_sel[-max_sel:]
            sel_idx = sel_idx[cls_sel]
            np.savetxt (walker_path + 'cls.sel.angle.0.out', sel_angles[cls_sel], fmt = '%.6f')
        elif len(sel_idx)==1:
            np.savetxt (walker_path + 'num_of_cluster.dat', [1], fmt = '%d')
        res_angles = np.loadtxt (walker_path + enhc_out_angle)
        res_angles = np.reshape (res_angles, [-1, cv_dim])
        res_angles = res_angles[sel_idx]
        np.savetxt (walker_path + 'cls.sel.out', sel_idx, fmt = '%d')
        np.savetxt (walker_path + 'cls.sel.angle.out', res_angles, fmt = '%.6f')
        res_confs = []
        for ii in sel_idx : 
            res_confs.append (walker_path + enhc_out_conf + ("conf%d.gro" % ii))    

        assert (len(res_confs) == res_angles.shape[0]), "number of enhc out conf does not match out angle"
        assert (len(sel_idx) == res_angles.shape[0]), "number of enhc out conf does not numb sel"
        nconf = len(res_confs)
        if nconf == 0 : 
            ret_list[walker_idx] = False
            continue

        sel_list=""
        for ii in range (nconf) : 
            if ii == 0 : sel_list = str(sel_idx[ii])
            else : sel_list += "," + str(sel_idx[ii])
        log_task ("selected %d confs, indexes: %s" % (nconf, sel_list))

        for ii in range (conf_start, nconf, conf_every) :
            # print (ii, sel_idx[ii])
            work_path = res_path + ((walker_format + ".%06d") % (walker_idx, sel_idx[ii])) + "/"
            os.makedirs(work_path)
            copy_file_list (mol_files, templ_mol_path, work_path)
            copy_file_list (res_files, templ_res_path, work_path)
            conf_file = walker_path + enhc_out_conf + ("conf%d.gro" % sel_idx[ii])
            if os.path.exists (work_path + "conf.gro") :
                os.remove (work_path + "conf.gro")
            conf_file = os.path.abspath(conf_file)
            tmp_cwd = os.getcwd()
            os.chdir(work_path)
            os.symlink (os.path.relpath(conf_file), "conf.gro")
            os.chdir(tmp_cwd)

        task_dirs = []
        task_args = []
        for ii in range (conf_start, nconf, conf_every) :
            dir_str = ((walker_format + ".%06d") % (walker_idx, sel_idx[ii]))
            arg_str = np.array2string(res_angles[ii], 
                                      formatter={'float_kind':lambda x: "%.6f" % x}).replace("[","").replace("]","").replace("\n"," ")
            task_dirs.append (dir_str)
            task_args.append (arg_str)
            log_task (task_dirs[-1] + ": " + task_args[-1])

        os.chdir (res_path)
        exec_hosts (MachineLocal, "./general_mkres.sh", 1, task_dirs, task_args)
        os.chdir (base_path)

        for ii in range (conf_start, nconf, conf_every) :
            work_path = res_path + ((walker_format + ".%06d") % (walker_idx, sel_idx[ii])) + "/"
            mol_conf_file = work_path + "grompp.mdp"
            make_grompp_res (mol_conf_file, nsteps, frame_freq)
            replace (work_path + res_plm, "STRIDE=[^ ]* ", "STRIDE=%d " % frame_freq)

    if any(ret_list) :
        return True
    else :
        return False

def run_res (iter_index,
             json_file, 
             exec_machine = MachineLocal) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    gmx_prep = jdata["gmx_prep"]
    gmx_run = jdata["gmx_run"]
    res_thread = jdata["res_thread"]
    gmx_run = gmx_run + (" -nt %d " % res_thread)
    gmx_run = gmx_run + " -plumed " + res_plm
    gmx_cont_run = gmx_run + " -cpi state.cpt "
    gmx_prep_log = "gmx_grompp.log"
    gmx_run_log = "gmx_mdrun.log"
    gmx_prep_cmd = cmd_append_log (gmx_prep, gmx_prep_log)
    gmx_run_cmd = cmd_append_log (gmx_run, gmx_run_log)
    gmx_cont_run_cmd = cmd_append_log (gmx_cont_run, gmx_run_log)
    res_group_size = jdata['res_group_size']
    batch_jobs = jdata['batch_jobs']
    batch_time_limit = jdata['batch_time_limit']
    batch_modules = jdata['batch_modules']
    batch_sources = jdata['batch_sources']
    
    iter_name = make_iter_name (iter_index)
    res_path = iter_name + "/" + res_name + "/"  
    base_path = os.getcwd() + "/"

    if not os.path.isdir (res_path) : 
        raise RuntimeError ("do not see any restrained simulation (%s)." % res_path)

    all_task_propose = glob.glob(res_path + "/[0-9]*[0-9]")
    if len(all_task_propose) == 0 :
        return
    all_task_propose.sort()
    if batch_jobs :
        all_task = all_task_propose
    else :
        all_task = []
        all_cont_task = []
        for ii in all_task_propose :
            if not os.path.isfile(os.path.join(ii, "confout.gro")) :
                if os.path.isfile(os.path.join(ii, "state.cpt")) :
                    all_cont_task.append(ii)
                else :
                    all_task.append(ii)

    if batch_jobs:
        exec_hosts (MachineLocal, gmx_prep_cmd, 1, all_task, None)
        exec_batch_group(gmx_run_cmd, res_thread, 1, all_task, task_args = None, group_size = res_group_size, time_limit = batch_time_limit, modules = batch_modules, sources = batch_sources)
    else :
        if len(all_task) == 1 :
            exec_hosts (MachineLocal, gmx_prep_cmd, 1, all_task, None)
            exec_hosts (MachineLocal, gmx_run_cmd, res_thread, all_task, None)
        elif len(all_task) > 1 :
            exec_hosts (MachineLocal, gmx_prep_cmd, 1, all_task, None)
            exec_hosts_batch (exec_machine, gmx_run_cmd, res_thread, all_task, None)
        if len(all_cont_task) == 1 :
            exec_hosts (MachineLocal, gmx_cont_run_cmd, res_thread, all_cont_task, None)
        elif len(all_cont_task) > 1 :
            exec_hosts_batch (exec_machine, gmx_cont_run_cmd, res_thread, all_cont_task, None)

def post_res (iter_index,
              json_file) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    res_cmpf_error = jdata["res_cmpf_error"]

    iter_name = make_iter_name (iter_index)
    res_path = iter_name + "/" + res_name + "/"  
    base_path = os.getcwd() + "/"
    
    all_task = glob.glob(res_path + "/[0-9]*[0-9]")
    if len(all_task) == 0 :
        np.savetxt (res_path + 'data.raw', [], fmt = "%.6e")        
        return
    all_task.sort()
    if res_cmpf_error :
        cmpf_cmd = "./cmpf.sh"
    else :
        cmpf_cmd = "./cmpf.py"
    cmpf_log = "cmpf.log"
    cmpf_cmd = cmd_append_log (cmpf_cmd, cmpf_log)

    centers = []
    force = []
    ndim = 0

    # run_node_tasks(max_thread, 1, all_task, cmpf_cmd)
    exec_hosts (MachineLocal, cmpf_cmd, 1, all_task, None)

    for work_path in all_task:
        os.chdir (work_path)
        this_centers = np.loadtxt ('centers.out')
        centers = np.append (centers, this_centers)
        this_force = np.loadtxt ('force.out')
        force = np.append (force, this_force)        
        ndim = this_force.size
        assert (ndim == this_centers.size), "center size is diff to force size in " + work_path
        os.chdir(base_path)        

    centers = np.reshape (centers, [-1, ndim])
    force = np.reshape (force, [-1, ndim])
    data = np.concatenate ((centers, force), axis = 1)
    np.savetxt (res_path + 'data.raw', data, fmt = "%.6e")

    norm_force = np.linalg.norm (force, axis = 1)
    log_task ("min|f| = %e  max|f| = %e  avg|f| = %e" % 
              (np.min(norm_force), np.max(norm_force), np.average(norm_force)))

def clean_res(iter_index):              
    iter_name = make_iter_name (iter_index)
    res_path = iter_name + "/" + res_name + "/"  
    base_path = os.getcwd() + "/"    

    all_task = glob.glob(res_path + "/[0-9]*[0-9]")
    all_task.sort()

    cleaned_files = ['cmpf*', '*log', 'general_mkres.sh', 'plm.res.out', 'state*.cpt', 'traj_comp.xtc', 'topol.tpr', 'tools', 'confout.gro', 'ener.edr', 'mdout.mdp', 'plumed.res.templ']
    cwd = os.getcwd()
    for ii in all_task :
        os.chdir(ii)
        clean_files(cleaned_files)
        os.chdir(cwd)    

def make_train (iter_index, 
                json_file) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    template_dir = jdata["template_dir"]    
    numb_model = jdata["numb_model"]
    res_iter = jdata["res_iter"]

    iter_name = make_iter_name (iter_index)
    train_path = iter_name + "/" + train_name + "/"  
    data_path = train_path + "data/"
    data_file = train_path + "data/data.raw"
    data_old_file = train_path + "data/data.old.raw"
    data_new_file = train_path + "data/data.new.raw"
    base_path = os.getcwd() + "/"
    templ_train_path = template_dir + "/" + train_name + "/"

    create_path (train_path) 
    os.makedirs(data_path)
    copy_file_list (train_files, templ_train_path, train_path)
    
    # collect data
    log_task ("collect data upto %d" % (iter_index))
    if iter_index == 0 :
        ii = 0
        this_raw = base_path + make_iter_name (ii) + "/" + res_name + "/data.raw"
        os.chdir(data_path)
        os.symlink (os.path.relpath(this_raw), os.path.basename(data_new_file))
        os.symlink (os.path.basename(data_new_file), os.path.basename(data_file))
        os.chdir (base_path)
        open (data_old_file, "w").close()
    else :
        prev_iter_index = iter_index - 1
        prev_data_file = base_path + make_iter_name(prev_iter_index) + "/" + train_name + "/data/data.raw"
        this_raw = base_path + make_iter_name (iter_index) + "/" + res_name + "/data.raw"
        os.chdir(data_path)
        os.symlink (os.path.relpath(prev_data_file), os.path.basename(data_old_file))
        os.symlink (os.path.relpath(this_raw), os.path.basename(data_new_file))
        os.chdir(base_path)
        with open(data_file, "wb") as fo :
            with open(data_old_file, "rb") as f0, open(data_new_file, "rb") as f1 :
                shutil.copyfileobj(f0, fo)
                shutil.copyfileobj(f1, fo)

    # create train dirs
    log_task ("create train dirs")
    for ii in range(numb_model) :
        work_path = train_path + ("%03d/" % ii) 
        old_model_path = work_path + "old_model/"
        create_path (work_path)
        os.chdir (work_path)
        os.symlink ("../data", "./data")
        os.chdir (base_path)
        if iter_index >= 1 :
            prev_iter_index = iter_index - 1
            prev_iter_name = make_iter_name (prev_iter_index)
            prev_train_path = prev_iter_name + "/" + train_name + "/"  
            prev_train_path = os.path.abspath(prev_train_path) + "/"
            prev_work_path = prev_train_path + ("%03d/" % ii)
            prev_model_files = glob.glob (prev_work_path + "model.ckpt.*")
            prev_model_files = prev_model_files + [prev_work_path + "checkpoint"]
            create_path (old_model_path)
            os.chdir(old_model_path)
            for ii in prev_model_files :
                os.symlink(os.path.relpath(ii), os.path.basename(ii))
                # shutil.copy (ii, old_model_path)
            os.chdir(base_path)
            for ii in prev_model_files :
                shutil.copy (ii, work_path)            

def run_train (iter_index, 
               json_file, 
               exec_machine = MachineLocal) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    numb_model = jdata["numb_model"]
    train_thread = jdata["train_thread"]
    res_iter = jdata["res_iter"]

    iter_name = make_iter_name (iter_index)
    train_path = iter_name + "/" + train_name + "/"  
    base_path = os.getcwd() + "/"

    # check if new data is empty
    new_data_file = os.path.join(train_path, 'data/data.new.raw')
    if os.stat(new_data_file).st_size == 0 :
        prev_iter_index = iter_index - 1
        prev_train_path = base_path + make_iter_name(prev_iter_index) + "/" + train_name + "/"
        prev_models = glob.glob(prev_train_path + "*.pb")
        for ii in prev_models :
            model_name = os.path.basename(ii)
            os.symlink(ii, os.path.join(train_path, model_name))
        return
    
    neurons = jdata["neurons"]
    batch_size = jdata["batch_size"]
    if iter_index < res_iter :
        numb_epoches = jdata["numb_epoches"]
        starter_lr = jdata["starter_lr"]
        decay_steps = jdata["decay_steps"]
        decay_rate = jdata["decay_rate"]    
        cmdl_args = ""
    else :
        numb_epoches = jdata["res_numb_epoches"]
        starter_lr = jdata["res_starter_lr"]
        decay_steps = jdata["res_decay_steps"]
        decay_rate = jdata["res_decay_rate"]            
        old_ratio = jdata["res_olddata_ratio"]
        cmdl_args = " --restart --use-mix --old-ratio %f " % old_ratio

    if jdata["resnet"] :
        cmdl_args += " --resnet "
    cmdl_args += " -n "
    for nn in neurons:
        cmdl_args += "%d " % nn
    cmdl_args += " -b " + str(batch_size)
    cmdl_args += " -e " + str(numb_epoches)
    cmdl_args += " -l " + str(starter_lr)
    cmdl_args += " --decay-steps " + str(decay_steps)
    cmdl_args += " --decay-rate " + str(decay_rate)

    train_cmd = "../main.py -t %d" % train_thread
    train_cmd += cmdl_args
    train_cmd = cmd_append_log (train_cmd, "train.log")
    freez_cmd = "../freeze.py -o graph.pb"
    freez_cmd = cmd_append_log (freez_cmd, "freeze.log")    
    task_dirs = [ ("%03d" % ii) for ii in range(numb_model) ]

    batch_jobs = jdata['batch_jobs']
    batch_time_limit = jdata['batch_time_limit']
    batch_modules = jdata['batch_modules']
    batch_sources = jdata['batch_sources']
    
    os.chdir (train_path)
    if batch_jobs:
        exec_batch(train_cmd, train_thread, 1, task_dirs, task_args = None, time_limit = batch_time_limit, modules = batch_modules, sources = batch_sources)
    else :
        if len(task_dirs) == 1 :
            exec_hosts (MachineLocal, train_cmd, train_thread, task_dirs, None)
        else :
            exec_hosts_batch (exec_machine, train_cmd, train_thread, task_dirs, None)

    exec_hosts (MachineLocal, freez_cmd, 1, task_dirs, None)
    for ii in range(numb_model) :
        os.symlink ("%03d/graph.pb" % ii, "graph.%03d.pb" % ii)
    os.chdir (base_path)

def clean_train(iter_index):              
    iter_name = make_iter_name (iter_index)
    train_path = iter_name + "/" + train_name + "/"  
    base_path = os.getcwd() + "/"    

    all_task = glob.glob(train_path + "/[0-9]*[0-9]")
    all_task.sort()

    # cleaned_files = ['checkpoint', 'model.ckpt*', 'freeze.log']
    cleaned_files = ['freeze.log']
    cwd = os.getcwd()
    for ii in all_task :
        os.chdir(ii)
        clean_files(cleaned_files)
        os.chdir(cwd)    
    
