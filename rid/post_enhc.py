
import os, json, glob, shutil
import numpy as np
from ridkit.lib.utils import make_iter_name, make_walker_name, cmd_append_log, set_resource, set_machine, log_task

from ridkit.lib.cal_cv_dim import cal_cv_dim
from ridkit.lib.cmpf import cmpf

from dpdispatcher.submission import Submission, Job, Task, Resources

enhc_out_conf="confs/"
enhc_out_angle="angle.rad.out"
enhc_out_plm="plm.out"
enhc_name="00.enhcMD"

res_name="01.resMD"
res_plm="plumed.res.dat"


def post_enhc (iter_index, 
               json_file,
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

    machine = set_machine(json_file)
    resources = set_resource(json_file, target="post")
    all_task_relpath = [os.path.relpath(ii, work_path) for ii in all_task]
    gmx_split_task = [ Task(command=gmx_split_cmd, task_work_path=ii, outlog='gmx_split.log', errlog='gmx_split.log') for ii in all_task_relpath ]
    gmx_split_submission = Submission(work_base=work_path, resources=resources, batch=machine, task_list=gmx_split_task)
    gmx_split_submission.run_submission()
    
    for ii in range(numb_walkers) :
        walker_path = work_path + make_walker_name(ii) + "/"
        angles = np.loadtxt (walker_path + enhc_out_plm)
        np.savetxt (walker_path + enhc_out_angle, angles[:,1:], fmt="%.6f")
    print("Post process of enhanced sampling finished.")