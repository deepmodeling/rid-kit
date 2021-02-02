#!/usr/bin/env python3

import os, re, time, logging, shutil
import numpy as np

iter_format = "%06d"
task_format = "%02d"
log_iter_head = "iter " + iter_format + " task " + task_format + ": "

def make_iter_name (iter_index) :
    return "iter." + (iter_format % iter_index)

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

def replace (file_name, pattern, subst) :
    file_handel = open (file_name, 'r')
    file_string = file_handel.read ()
    file_handel.close ()
    file_string = ( re.sub (pattern, subst, file_string) )
    file_handel = open (file_name, 'w')
    file_handel.write (file_string)
    file_handel.close ()

def make_grompp_res (gro_file, nsteps, frame_freq) :
    replace (gro_file, "nsteps.*=.*", "nsteps = %d" % nsteps)
    replace (gro_file, "nstxout.*=.*", "nstxout = %d" % 0)
    replace (gro_file, "nstvout.*=.*", "nstvout = %d" % 0)
    replace (gro_file, "nstfout.*=.*", "nstfout = %d" % 0)
    replace (gro_file, "nstenergy.*=.*", "nstenergy = %d" % 0)

def make_grompp_res_equi (gro_file, nsteps, frame_freq, dt) :
    replace (gro_file, "nsteps.*=.*", "nsteps = %d" % nsteps)
    replace (gro_file, "dt.*=.*", "dt = %f" % dt)
    replace (gro_file, "nstxout.*=.*", "nstxout = %d" % 0)
    replace (gro_file, "nstvout.*=.*", "nstvout = %d" % 0)
    replace (gro_file, "nstfout.*=.*", "nstfout = %d" % 0)
    replace (gro_file, "nstenergy.*=.*", "nstenergy = %d" % 0)

def copy_file_list (file_list, from_path, to_path) :
    for jj in file_list : 
        if os.path.isfile(from_path + jj) :
            shutil.copy (from_path + jj, to_path)
        elif os.path.isdir(from_path + jj) :
            shutil.copytree (from_path + jj, to_path + jj)

def create_path (path) :
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

def list_to_arg (mylist) :
    return np.array2string(mylist,
                           formatter={'float_kind':lambda x: "%.6f" % x}).replace("[","").replace("]","").replace("\n"," ")    

def exec_hosts (machine_env,
                cmd, 
                work_thread,
                task_dirs,
                task_args = None) :
    ntasks = len(task_dirs)
    if task_args != None :
        assert ntasks == len(task_args) or len(task_args) == 1
        if len(task_args) == 1 :
            tmp_arg = task_args[0]
            task_args = [tmp_arg for ii in range(ntasks)]
    else :
        task_args = ["" for ii in range(ntasks)]
    assert ntasks == len(task_args)    

    host_list = machine_env.get_node_list()
    nnode = len(host_list)
    ntpnode = machine_env.get_core_per_node()
    nsource = ntpnode * nnode
    numb_jobs = (nsource // work_thread)
    task_chunks = [task_dirs[i:i + numb_jobs] for i in range(0, ntasks, numb_jobs)]
    args_chunks = [task_args[i:i + numb_jobs] for i in range(0, ntasks, numb_jobs)]
    nbatch = len(task_chunks)
    assert nbatch == len(args_chunks)

    base_path = os.getcwd() + "/"
    for ii in range(nbatch) :
        task_batch = task_chunks[ii]
        args_batch = args_chunks[ii]
        ps = []
        for jj in range(len(task_batch)) :
            work_path = task_batch[jj]
            work_args = args_batch[jj]
            host = host_list[jj % nnode]
            os.chdir(work_path)
            work_name = os.path.basename (os.getcwd())
            os.chdir(base_path)
            logging.info(("%s %03d %s: %s %s" % (host, ii, work_name, cmd, work_args)))
            ps.append(machine_env.exec_cmd(host, cmd, work_path, work_args))
        while True :
            if not(any(p.wait() for p in ps)) :
                break
            time.sleep(1)    


