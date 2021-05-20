import os, json, glob
from dpdispatcher.lazy_local_context import LazyLocalContext
from dpdispatcher.submission import Submission, Job, Task, Resources

from rid.lib.utils import make_iter_name, make_walker_name, cmd_append_log, set_resource, set_machine
from rid.lib.cal_cv_dim import cal_cv_dim

res_plm="plumed.res.dat"
res_name="01.resMD"

def run_res (iter_index,
             json_file,
             base_dir="./") :
    json_file = os.path.abspath(json_file)
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    fp.close()
    gmx_prep = jdata["gmx_prep"]
    gmx_run = jdata["gmx_run"]
    res_thread = jdata["res_thread"]
    gmx_run = gmx_run + (" -nt %d" % res_thread)
    gmx_run = gmx_run + " -plumed " + res_plm
    # gmx_cont_run = gmx_run + " -cpi state.cpt"
    gmx_prep_log = "gmx_grompp.log"
    gmx_run_log = "gmx_mdrun.log"
    gmx_prep_cmd = cmd_append_log (gmx_prep, gmx_prep_log)
    gmx_run_cmd = cmd_append_log (gmx_run, gmx_run_log)
    # gmx_cont_run_cmd = cmd_append_log (gmx_cont_run, gmx_run_log)
    
    base_dir = os.path.abspath(base_dir) + "/"
    iter_name = make_iter_name (iter_index)
    res_path = base_dir + iter_name + "/" + res_name + "/"  

    if not os.path.isdir (res_path) : 
        raise RuntimeError ("do not see any restrained simulation (%s)." % res_path)
    
    all_task = list(filter(lambda x:os.path.isdir(x),  glob.glob(res_path + "/[0-9]*[0-9]")))
    print('run_res:all_task_propose:', all_task)
    print('run_res:gmx_prep_cmd:', gmx_prep_cmd)
    print('run_res:gmx_run_cmd:', gmx_run_cmd)
    # print('run_res:gmx_cont_run_cmd:', gmx_cont_run_cmd)

    if len(all_task) == 0:
        return None
    all_task.sort()
    all_task_basedir = [os.path.relpath(ii, res_path) for ii in all_task]

    res_resources = set_resource(json_file, target="res")
    machine = set_machine(json_file)

    gmx_prep_task = [ Task(command=gmx_prep_cmd, task_work_path=ii, outlog='gmx_grompp.log', errlog='gmx_grompp.log') for ii in all_task_basedir ]
    gmx_prep_submission = Submission(work_base=res_path, resources=res_resources, batch=machine, task_list=gmx_prep_task)
    gmx_prep_submission.run_submission()

    gmx_run_task =  [ Task(command=gmx_run_cmd, task_work_path=ii, outlog='gmx_mdrun.log', errlog='gmx_mdrun.log') for ii in all_task_basedir ]
    gmx_run_submission = Submission(work_base=res_path, resources=res_resources, batch=machine, task_list=gmx_run_task)
    gmx_run_submission.run_submission()

