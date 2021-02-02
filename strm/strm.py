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

exec_machine = MachineLocal

template_path = os.path.abspath("./template")

res_name="00.resMD"
res_files=["cmpf.sh", "cmpf.py", "mkres.sh", "plumed.res.templ", "tools"]
res_plm = "plumed.res.dat"

mol_name="mol"
mol_files=["grompp.mdp", "topol.top", "angle.ndx", "angle.sh"]
conf_angle_files = ["angle.ndx", "conf.ang.sh"]
string_equi_conf_files = ["grompp.mdp", "topol.top"]

def cvt_conf_angle (conf_file) :
    base_path = os.getcwd() + "/"
    mol_path = template_path + "/" + mol_name + "/"
    conf_file_path = os.path.abspath (conf_file)

    conf_angle_path = base_path + "conf_angle_tmp" + "/"
    libr.create_path (conf_angle_path)
    libr.copy_file_list (conf_angle_files, mol_path, conf_angle_path)
    
    os.chdir (conf_angle_path)
    sp.check_call ("./conf.ang.sh " + conf_file_path + " > /dev/null", shell = True)
    ang = np.loadtxt ("angle.rad.out")
    os.chdir (base_path)
    shutil.rmtree (conf_angle_path)
    
    return ang

def init_linear_string (node_start, node_end, numb_node):
    """ init a linear string between the start and end nodes. """
    dim = np.size(node_start)
    if dim != np.size(node_end):
        raise NameError ('The dimention of starting and ending nodes of the string should match!')
    string = np.zeros ((dim, numb_node))
    for ii in range (dim):
        string [ii] = np.linspace (node_start[ii], node_end[ii], numb_node)
    return string.T

def make_string (iter_index,
                 string) :
    iter_name = libr.make_iter_name (iter_index)
    base_path = os.getcwd() + "/"
    string_path = base_path + iter_name + "/"
    libr.create_path (string_path)

    os.chdir (string_path)
    np.savetxt ("string.out", string)
    os.chdir (base_path)

def equi_0_serial (json_file) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    equi_nsteps = jdata["res_equi_nsteps"]
    equi_dt = jdata["res_equi_dt"]
    init_start_conf = jdata["init_start_conf"]
    start_conf_name = os.path.basename(init_start_conf)
    gmx_prep = jdata["gmx_prep"]
    gmx_run = jdata["gmx_run"]
    iter_index = 0
    frame_freq = 0
    gmx_prep_log = "gmx_grompp.log"
    gmx_run_log = "gmx_mdrun.log"

    base_path = os.getcwd() + "/"
    iter_name = libr.make_iter_name (iter_index)
    string_path = base_path + iter_name + "/"
    res_path = string_path + res_name + ".equi" + "/"
    templ_mol_path = template_path + "/" + mol_name + "/"
    templ_res_path = template_path + "/" + res_name + "/"

    libr.create_path (res_path) 
    
    string = np.loadtxt (string_path + "string.out")
    numb_node = string.shape[0]
        
    # case of node 0
    node_name = "%06d" % 0
    work_path = res_path + node_name + "/"
    libr.create_path(work_path)
    libr.copy_file_list(res_files,               templ_res_path, work_path)    
    libr.copy_file_list(string_equi_conf_files,  templ_mol_path, work_path)
    libr.copy_file_list([init_start_conf],            base_path, work_path)
    os.chdir(work_path)
    if not os.path.isfile("conf.gro") :
        os.symlink(start_conf_name, "conf.gro")
    os.chdir(base_path)
    ## equi simul
    os.chdir (work_path)
    arg_str = libr.list_to_arg(string[0])
    sp.check_call("./mkres.sh " + arg_str, shell = True)
    libr.log_task ("./mkres.sh " + arg_str)
    libr.make_grompp_res_equi("grompp.mdp", equi_nsteps, frame_freq, equi_dt)
    libr.replace (work_path + res_plm, "STRIDE=[^ ]* ", "STRIDE=%d " % frame_freq)
    gmx_prep_cmd = gmx_prep 
    gmx_run_cmd = gmx_run + " -nt 1 " + " -plumed " + res_plm
    gmx_prep_cmd = libr.cmd_append_log (gmx_prep_cmd, gmx_prep_log)
    gmx_run_cmd = libr.cmd_append_log (gmx_run_cmd, gmx_prep_log)
    libr.log_task(gmx_prep_cmd)
    sp.check_call(gmx_prep_cmd, shell = True)
    libr.log_task(gmx_run_cmd)
    sp.check_call(gmx_run_cmd, shell = True)
    os.chdir(base_path)
    
    for ii in range (1, numb_node) :
        prev_node_name = "%06d" % (ii-1)
        node_name = "%06d" % ii
        work_path = res_path + node_name + "/"
        libr.create_path(work_path)
        libr.copy_file_list(res_files,               templ_res_path, work_path)    
        libr.copy_file_list(string_equi_conf_files,  templ_mol_path, work_path)
        os.symlink(work_path + "../" + prev_node_name + "/confout.gro", work_path + "/conf.gro")
        ## equi simul
        os.chdir (work_path)
        arg_str = libr.list_to_arg(string[ii])
        libr.log_task("./mkres.sh " + arg_str)
        sp.check_call("./mkres.sh " + arg_str, shell = True)
        libr.make_grompp_res_equi("grompp.mdp", equi_nsteps, frame_freq, equi_dt)
        libr.replace(work_path + res_plm, "STRIDE=[^ ]* ", "STRIDE=%d " % frame_freq)
        gmx_prep_cmd = gmx_prep
        gmx_run_cmd = gmx_run + " -nt 1 " + " -plumed " + res_plm
        gmx_prep_cmd = libr.cmd_append_log(gmx_prep_cmd, gmx_prep_log)
        gmx_run_cmd = libr.cmd_append_log(gmx_run_cmd, gmx_prep_log)
        libr.log_task(gmx_prep_cmd)
        sp.check_call(gmx_prep_cmd, shell = True)
        libr.log_task(gmx_run_cmd)
        sp.check_call(gmx_run_cmd, shell = True)
        os.chdir(base_path)        

def equi_old_string_parallel (json_file) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    equi_nsteps = jdata["res_equi_nsteps"]
    equi_dt = jdata["res_equi_dt"]
    old_string_path = jdata["old_string_path"]
    assert os.path.isdir(old_string_path)
    iter_index = 0
    frame_freq = 0

    old_string_path = os.path.abspath(old_string_path) + "/"
    base_path = os.getcwd() + "/"
    iter_name = libr.make_iter_name (iter_index)
    new_string_path = base_path + iter_name + "/"
    res_path = new_string_path + res_name + ".equi" + "/"
    templ_mol_path = template_path + "/" + mol_name + "/"
    templ_res_path = template_path + "/" + res_name + "/"

    old_string = np.loadtxt(old_string_path + "string.out")
    new_string = np.loadtxt(new_string_path + "string.out")
    old_numb_nodes = old_string.shape[0]
    new_numb_nodes = new_string.shape[0]
    
    for ii in range(new_numb_nodes) :
        # create working dir
        node_name = "%06d" % ii
        work_path = res_path + node_name + "/"
        libr.create_path(work_path)
        libr.copy_file_list(res_files,               templ_res_path, work_path)    
        libr.copy_file_list(string_equi_conf_files,  templ_mol_path, work_path)
        # find the minimal dist old node
        min_indx = 0
        min_dist = np.linalg.norm(new_string[ii] - old_string[min_indx])
        for jj in range(1, old_numb_nodes):
            dist = np.linalg.norm(new_string[ii] - old_string[jj])
            if dist < min_dist :
                min_indx = jj
                min_dist = dist
        init_node_path = old_string_path + res_name + "/" + ("%06d"%min_indx) + "/"
        # link confs
        os.symlink(init_node_path + "/confout.gro", work_path + "/conf.gro")
        # setup working dir
        os.chdir (work_path)
        arg_str = libr.list_to_arg(new_string[ii])
        libr.log_task("./mkres.sh " + arg_str)
        sp.check_call("./mkres.sh " + arg_str, shell = True)
        libr.make_grompp_res_equi("grompp.mdp", equi_nsteps, frame_freq, equi_dt)
        libr.replace(work_path + res_plm, "STRIDE=[^ ]* ", "STRIDE=%d " % frame_freq)
        os.chdir(base_path)        

    run_res(0, json_file, ".equi")

def make_res_x (iter_index, 
                json_file) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    nsteps = jdata["res_nsteps"]
    gmx_prep = jdata["gmx_prep"]
    gmx_run = jdata["gmx_run"]
    frame_freq = jdata["res_frame_freq"]

    base_path = os.getcwd() + "/"
    string_path      = base_path + libr.make_iter_name(iter_index)   + "/"
    string_equi_path = base_path + libr.make_iter_name(iter_index-1) + "/"
    res_path = string_path + res_name + "/"
    if iter_index == 0:
        res_equi_path = string_path + res_name + ".equi" + "/"
    else :
        res_equi_path = string_equi_path + res_name + "/"
    templ_mol_path = template_path + "/" + mol_name + "/"
    templ_res_path = template_path + "/" + res_name + "/"

    libr.create_path (res_path) 
    
    string = np.loadtxt (string_path + "string.out")
    numb_node = string.shape[0]

    task_dirs = []
    task_args = []
    for ii in range (0, numb_node) :
        node_name = "%06d" % ii
        work_path = res_path + node_name + "/"
        equi_path = res_equi_path + node_name + "/"
        libr.create_path(work_path)
        libr.copy_file_list(res_files, templ_res_path, work_path)    
        libr.copy_file_list(mol_files, templ_mol_path, work_path)
        os.symlink(equi_path + "confout.gro", work_path + "conf.gro")
        ## simul
        arg_str = libr.list_to_arg(string[ii])
        # libr.log_task ("./mkres.sh " + arg_str)
        # sp.check_call("./mkres.sh " + arg_str, shell = True)
        task_dirs.append(work_path)
        task_args.append(arg_str)

    global exec_machine
    libr.exec_hosts (exec_machine, "./mkres.sh", 1, task_dirs, task_args)

    for ii in range (0, numb_node) :
        node_name = "%06d" % ii
        work_path = res_path + node_name + "/"
        os.chdir (work_path)
        libr.make_grompp_res("grompp.mdp", nsteps, frame_freq)
        libr.replace (work_path + res_plm, "STRIDE=[^ ]* ", "STRIDE=%d " % frame_freq)    
        os.chdir(base_path)

def make_res (iter_index, 
              json_file): 
    make_res_x (iter_index, json_file)

def run_res (iter_index,
             json_file, 
             res_name_suffix = "") :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    gmx_prep = jdata["gmx_prep"]
    gmx_run = jdata["gmx_run"]
    res_thread = jdata["res_thread"]
    gmx_run = gmx_run + (" -nt %d " % res_thread)
    gmx_run = gmx_run + " -plumed " + res_plm
    gmx_prep_log = "gmx_grompp.log"
    gmx_run_log = "gmx_mdrun.log"
    gmx_prep_cmd = libr.cmd_append_log (gmx_prep, gmx_prep_log)
    gmx_run_cmd = libr.cmd_append_log (gmx_run, gmx_run_log)
    
    iter_name = libr.make_iter_name (iter_index)
    res_path = iter_name + "/" + res_name + res_name_suffix + "/"  
    base_path = os.getcwd() + "/"

    if not os.path.isdir (res_path) : 
        raise RuntimeError ("do not see any restrained simulation (%s)." % res_path)

    all_task = glob.glob(res_path + "/[0-9]*[0-9]")
    all_task.sort()

    # run_node_tasks(max_thread, 1, all_task, gmx_prep_cmd)
    # run_node_tasks(max_thread, res_thread, all_task, gmx_run_cmd)
    global exec_machine
    libr.exec_hosts (exec_machine, gmx_prep_cmd, 1, all_task, None)
    libr.exec_hosts (exec_machine, gmx_run_cmd, res_thread, all_task, None)

def post_res (iter_index,
              json_file) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    res_cmpf_error = jdata["res_cmpf_error"]

    iter_name = libr.make_iter_name (iter_index)
    res_path = iter_name + "/" + res_name + "/"  
    base_path = os.getcwd() + "/"
    string_path = base_path + iter_name + "/"
    if not os.path.isdir (res_path) : 
        raise RuntimeError ("do not see any restrained simulation (%s)." % res_path)
    
    all_task = glob.glob(res_path + "/[0-9]*[0-9]")
    all_task.sort()
    if res_cmpf_error :
        cmpf_cmd = "./cmpf.sh"
    else :
        cmpf_cmd = "./cmpf.py"
    cmpf_log = "cmpf.log"
    cmpf_cmd = libr.cmd_append_log (cmpf_cmd, cmpf_log)

    global exec_machine
    libr.exec_hosts (exec_machine, cmpf_cmd, 1, all_task, None)

    centers = []
    force = []
    ndim = 0
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
    # added this one...
    np.savetxt (string_path + 'force.out', force, fmt = "%.6e")

    norm_force = np.linalg.norm (force, axis = 1)
    libr.log_task ("min|f| = %e  max|f| = %e  avg|f| = %e" % 
                   (np.min(norm_force), np.max(norm_force), np.average(norm_force)))    

def compute_force (iter_index,
                   string_) :
    base_path = os.getcwd() + "/"
    iter_name = libr.make_iter_name (iter_index)
    string_path = base_path + iter_name + "/"    
    string = np.loadtxt (string_path + "string.out")
    norm_diff = np.linalg.norm (string_ - string)
    if (norm_diff > 1e-5) :
        raise RuntimeError ("inconsistent string")
    return np.loadtxt (string_path + "force.out")

def update_string (iter_index, 
                   json_file) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    dt = jdata["str_dt"]
    weighting = [[0, 1], [1, 1]]
    key = "string_weight"
    if key in jdata.keys() :
       weighting = jdata[key] 

    base_path = os.getcwd() + "/"
    iter_name = libr.make_iter_name (iter_index)
    string_path = base_path + iter_name + "/"
    string = np.loadtxt (string_path + "string.out")
    numb_node = string.shape[0]
    
    old_string = libs.resample_string (string, numb_node, weighting)
    new_string = libs.update_string_Euler (compute_force, dt, iter_index, string)
    new_string = libs.resample_string (new_string, numb_node, weighting)
    
    make_string (iter_index + 1, new_string)
    
    return np.linalg.norm (old_string - new_string)

def init_string_confs_linear (json_file) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    init_start_conf = jdata["init_start_conf"]
    init_end_conf = jdata["init_end_conf"]
    numb_discr = jdata["numb_discr"]
    numb_nodes = numb_discr + 1

    startn = (cvt_conf_angle (init_start_conf))
    endn =   (cvt_conf_angle (init_end_conf))

    ndim = len(startn)
    for ii in range(ndim) :
        diff = startn[ii] - endn[ii]
        if  diff < - np.pi:
            startn[ii] += 2. * np.pi
        elif diff >= np.pi:
            startn[ii] -= 2. * np.pi

    string = init_linear_string (startn, endn, numb_nodes)
    ii = 0
    make_string(ii, string)

def init_string_resample_old (json_file) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)    
    old_string_path = jdata["old_string_path"]
    assert os.path.isdir(old_string_path)
    old_string_path = os.path.abspath(old_string_path) + "/"
    numb_discr = jdata["numb_discr"]
    numb_nodes = numb_discr + 1
    weighting = [[0, 1], [1, 1]]
    key = "string_weight"
    if key in jdata.keys() :
       weighting = jdata[key] 

    old_string = np.loadtxt(old_string_path + "string.out")
    string = libs.resample_string(old_string, numb_nodes, weighting)
    ii = 0
    make_string(ii, string)

def run_iter (json_file) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    numb_iter = jdata["numb_iter"]
    tol = jdata["tol_diff"]
    init_method = jdata["init_method"]
    numb_node = jdata["numb_discr"] + 1
    record = "record.strm"
    numb_task = 5
    max_tasks = 1000000

    iter_rec = [0, -1]
    if os.path.isfile (record) :
        with open (record) as frec :
            for line in frec : 
                iter_rec = [int(x) for x in line.split()]
        logging.info ("continue from iter %03d task %02d" % (iter_rec[0], iter_rec[1]))

    for ii in range (numb_iter) :
        for jj in range (numb_task) :
            if ii * max_tasks + jj <= iter_rec[0] * max_tasks + iter_rec[1] : 
                continue
            if jj == 0 :
                if ii == 0:         
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
            elif jj == 1:
                libr.log_iter ("make_res", ii, jj)
                make_res (ii, json_file)
            elif jj == 2 :
                libr.log_iter ("run_res", ii, jj)
                run_res (ii, json_file)
            elif jj == 3 :
                libr.log_iter ("post_res", ii, jj)
                post_res (ii, json_file)
            elif jj == 4 :
                diff = update_string (ii, json_file)
                diff /= np.sqrt(numb_node)
                libr.log_iter (("update_string with diff %e" % diff), ii, jj)
                if diff < tol :
                    libr.log_iter ("string converged", ii, jj)
                    return
            else :
                raise RuntimeError ("unknow task %d, something wrong" % jj)
            
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
