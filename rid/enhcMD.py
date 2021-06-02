import os, json, glob, shutil
import numpy as np
from rid.lib.utils import create_path
from rid.lib.utils import replace
from rid.lib.utils import make_iter_name
from rid.lib.utils import make_walker_name
from rid.lib.utils import checkfile
from rid.lib.utils import cmd_append_log, log_task

from rid.lib.gen.gen_mdp import make_grompp
from rid.lib.gen.gen_plumed import make_plumed
from rid.lib.gen.gen_plumed import conf_enhc_plumed
from rid.lib.cal_cv_dim import cal_cv_dim
from rid.lib.cmpf import cmpf

from rid.lib.machine import set_resource, set_machine
from dpdispatcher.submission import Submission, Job, Task

enhc_name="00.enhcMD"
enhc_plm="plumed.dat"
enhc_bf_plm="plumed.bf.dat"
enhc_out_plm="plm.out"
enhc_out_conf="confs/"
enhc_out_angle="angle.rad.out"

res_name="01.resMD"
res_plm="plumed.res.dat"

def adjust_lvl(prev_enhc_path, num_of_cluster_threshhold, jdata):
    # adaptive trust level
    num_of_cluster=np.loadtxt(prev_enhc_path+'num_of_cluster.dat')
    pre_trust_lvl1=np.loadtxt(prev_enhc_path+'trust_lvl1.dat')
    if num_of_cluster < num_of_cluster_threshhold:
        enhc_trust_lvl_1 = pre_trust_lvl1 * 1.5
        enhc_trust_lvl_2 = enhc_trust_lvl_1+1
    else:
        enhc_trust_lvl_1 = jdata["bias_trust_lvl_1"]
        enhc_trust_lvl_2 = enhc_trust_lvl_1+1
    if enhc_trust_lvl_1>jdata["bias_trust_lvl_1"]*8:
        enhc_trust_lvl_1 = jdata["bias_trust_lvl_1"]
        enhc_trust_lvl_2 = enhc_trust_lvl_1+1
    return enhc_trust_lvl_1, enhc_trust_lvl_2

def prep_graph(graph_files, walker_path):
    # copy graph files
    cwd = os.getcwd()
    for ii in graph_files :
        file_name = os.path.basename(ii)
        abs_path = os.path.abspath(ii)   
        rel_path = os.path.relpath(abs_path, os.path.abspath(walker_path))     
        checkfile(walker_path + file_name)
        os.chdir(walker_path)
        os.symlink(rel_path, walker_path + file_name)    
        os.chdir(cwd)

def get_graph_list(graph_files):
    graph_list=""
    counter=0
    for ii in graph_files :
        file_name = os.path.basename(ii)
        if counter == 0 :
            graph_list="%s" % file_name
        else :
            graph_list="%s,%s" % (graph_list, file_name)
        counter = counter + 1
    return graph_list

def make_enhc (iter_index, 
               json_file, 
               graph_files,
               mol_dir,
               cv_file,
               base_dir='./') :
    base_dir = os.path.abspath(base_dir) + "/"
    json_file = os.path.abspath(json_file)
    cv_file = os.path.abspath(cv_file)
    graph_files.sort()
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    fp.close()
    numb_walkers = jdata["numb_walkers"]
    enhc_trust_lvl_1 = jdata["bias_trust_lvl_1"]
    enhc_trust_lvl_2 = jdata["bias_trust_lvl_2"]
    nsteps = jdata["bias_nsteps"]
    frame_freq = jdata["bias_frame_freq"]
    num_of_cluster_threshhold = jdata["num_of_cluster_threshhold"]
    dt = jdata["bias_dt"]
    temperature = jdata["bias_temperature"]
    
    iter_name = make_iter_name(iter_index)
    work_path = base_dir + iter_name + "/" + enhc_name + "/"  
    mol_path = os.path.abspath(mol_dir) + "/"
    
    conf_list = glob.glob(mol_path + "*gro")
    conf_list.sort()
    assert (len(conf_list) >= numb_walkers), "not enough conf files in mol dir %s" % mol_path

    create_path(work_path)

    mol_files=["topol.top"]
    for walker_idx in range (numb_walkers) :
        walker_path = work_path + make_walker_name(walker_idx) + "/"
        create_path (walker_path)
       
        make_grompp(walker_path + "grompp.mdp", "bias", nsteps, frame_freq, temperature=temperature, dt=dt, define='')
        # make_grompp(walker_path + "grompp_restraint.mdp", "res", nsteps, frame_freq, temperature=temperature, dt=dt, define='-DPOSRE')
        
        for ii in mol_files:
            checkfile(walker_path + ii)
            shutil.copy(mol_path + ii, walker_path)

        # copy conf file
        conf_file = conf_list[walker_idx]
        checkfile(walker_path + "conf.gro")
        shutil.copy(conf_file, walker_path + "conf.gro")
        checkfile(walker_path + "conf_init.gro")
        shutil.copy(conf_file, walker_path + "conf_init.gro")

        # if have prev confout.gro, use as init conf
        if iter_index > 0 :
            prev_enhc_path = base_dir + make_iter_name(iter_index-1) + "/" + enhc_name + "/" + make_walker_name(walker_idx) + "/"
            prev_enhc_path = os.path.abspath(prev_enhc_path) + "/"
            if os.path.isfile (prev_enhc_path + "confout.gro") :
                os.remove (walker_path + "conf.gro")
                rel_prev_enhc_path = os.path.relpath(prev_enhc_path + "confout.gro", walker_path)
                os.symlink (rel_prev_enhc_path, walker_path + "conf.gro")
            else :
                raise RuntimeError("cannot find prev output conf file  " + prev_enhc_path + 'confout.gro')
            log_task ("use conf of iter " + make_iter_name(iter_index-1) + " walker " + make_walker_name(walker_idx) )

            enhc_trust_lvl_1, enhc_trust_lvl_2 = adjust_lvl(prev_enhc_path, num_of_cluster_threshhold, jdata)
            
        np.savetxt(walker_path+'trust_lvl1.dat', [enhc_trust_lvl_1] , fmt = '%.6f')

        make_plumed(walker_path, "dpbias", conf_file, cv_file)
        make_plumed(walker_path, "bf", conf_file, cv_file)

        prep_graph(graph_files, walker_path)
        # config plumed
        graph_list = get_graph_list(graph_files)
        conf_enhc_plumed(walker_path + enhc_plm, "enhc", graph_list, enhc_trust_lvl_1=enhc_trust_lvl_1, enhc_trust_lvl_2=enhc_trust_lvl_2, frame_freq=frame_freq, enhc_out_plm=enhc_out_plm)
        conf_enhc_plumed(walker_path + enhc_bf_plm, "bf", graph_list, frame_freq=frame_freq, enhc_out_plm=enhc_out_plm)
        
        if len(graph_list) == 0 :
            log_task ("brute force MD without NN acc")
        else :
            log_task ("use NN model(s): " + graph_list)
            log_task ("set trust l1 and l2: %f %f" % (enhc_trust_lvl_1, enhc_trust_lvl_2))
    print("Enhanced sampling has prepared.")


def run_enhc (iter_index,
              json_file,
              machine_json,
              base_dir='./') :
    json_file = os.path.abspath(json_file)
    base_dir = os.path.abspath(base_dir) + "/"
    iter_name = make_iter_name (iter_index)
    work_path = base_dir + iter_name + "/" + enhc_name + "/"  

    fp = open (json_file, 'r')
    jdata = json.load (fp)
    fp.close()
    gmx_prep = jdata["gmx_prep"]
    gmx_run = jdata["gmx_run"]
    enhc_thread = jdata["enhc_thread"]
    gmx_run = gmx_run + (" -nt %d" % enhc_thread)
    gmx_prep_log = "gmx_grompp.log"
    gmx_run_log = "gmx_mdrun.log"
    # assuming at least one walker
    graph_files = glob.glob (work_path + (make_walker_name(0)) + "/*.pb")
    if len (graph_files) != 0 :
        gmx_run = gmx_run + " -plumed " + enhc_plm  
    else :
        gmx_run = gmx_run + " -plumed " + enhc_bf_plm
    gmx_prep_cmd = cmd_append_log (gmx_prep, gmx_prep_log)
    gmx_run_cmd = cmd_append_log (gmx_run, gmx_run_log)
    numb_walkers = jdata["numb_walkers"]

    all_task = list(filter(lambda x:os.path.isdir(x),  glob.glob(work_path + "/[0-9]*[0-9]")))
    all_task.sort()

    all_task_basedir = [os.path.relpath(ii, work_path) for ii in all_task]
    print('run_enhc:work_path', work_path)
    print('run_enhc:gmx_prep_cmd:', gmx_prep_cmd)
    print('run_enhc:gmx_run_cmd:', gmx_run_cmd)
    print('run_enhc:all_task:', all_task)
    print('run_enhc:all_task_basedir:', all_task_basedir)
    
    machine = set_machine(machine_json, target="enhcMD")
    resources = set_resource(machine_json, target="enhcMD")

    gmx_prep_task = [ Task(command=gmx_prep_cmd, task_work_path=ii, outlog='gmx_grompp.log', errlog='gmx_grompp.log') for ii in all_task_basedir ]
    gmx_prep_submission = Submission(work_base=work_path, machine=machine, resources=resources, task_list=gmx_prep_task)

    gmx_prep_submission.run_submission()
    
    gmx_run_task =  [ Task(command=gmx_run_cmd, task_work_path=ii, outlog='gmx_mdrun.log', errlog='gmx_mdrun.log') for ii in all_task_basedir ]
    gmx_run_submission = Submission(work_base=work_path, machine=machine, resources=resources, task_list=gmx_run_task)
    gmx_run_submission.run_submission()


def post_enhc (iter_index, 
               json_file,
               machine_json,
               base_dir="./") :
    base_dir = os.path.abspath(base_dir) + "/"
    iter_name = make_iter_name (iter_index)
    work_path = base_dir + iter_name + "/" + enhc_name + "/" 
    json_file = os.path.abspath(json_file) 
    json_file = os.path.abspath(json_file)
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    fp.close()
    gmx_split = jdata["gmx_split_traj"]
    gmx_split_log = "gmx_split.log"
    gmx_split_cmd = cmd_append_log (gmx_split, gmx_split_log)
    
    all_task = list(filter(lambda x:os.path.isdir(x),  glob.glob(work_path + "/[0-9]*[0-9]")))
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

    print('rid.py:post_enhc:gmx_split_cmd', gmx_split_cmd)
    print('rid.py:post_enhc:work path', work_path)

    machine = set_machine(machine_json, target="post_enhc")
    resources = set_resource(machine_json, target="post_enhc")
    all_task_relpath = [os.path.relpath(ii, work_path) for ii in all_task]
    gmx_split_task = [ Task(command=gmx_split_cmd, task_work_path=ii, outlog='gmx_split.log', errlog='gmx_split.log') for ii in all_task_relpath ]
    gmx_split_submission = Submission(work_base=work_path, resources=resources, machine=machine, task_list=gmx_split_task)
    gmx_split_submission.run_submission()
    
    for ii in range(numb_walkers) :
        walker_path = work_path + make_walker_name(ii) + "/"
        angles = np.loadtxt (walker_path + enhc_out_plm)
        np.savetxt (walker_path + enhc_out_angle, angles[:,1:], fmt="%.6f")
    print("Post process of enhanced sampling finished.")
