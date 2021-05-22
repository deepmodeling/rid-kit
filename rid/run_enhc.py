import json, os, glob
from dpdispatcher.lazy_local_context import LazyLocalContext
from dpdispatcher.submission import Submission, Job, Task, Resources
from ridkit.lib.utils import make_iter_name, make_walker_name, cmd_append_log, set_resource, set_machine

enhc_name="00.enhcMD"
enhc_plm="plumed.dat"
enhc_bf_plm="plumed.bf.dat"

def run_enhc (iter_index,
              json_file,
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
    
    machine = set_machine(json_file)
    resources = set_resource(json_file, target="enhc")

    gmx_prep_task = [ Task(command=gmx_prep_cmd, task_work_path=ii, outlog='gmx_grompp.log', errlog='gmx_grompp.log') for ii in all_task_basedir ]
    gmx_prep_submission = Submission(work_base=work_path, resources=resources, batch=machine, task_list=gmx_prep_task)

    gmx_prep_submission.run_submission()
    
    gmx_run_task =  [ Task(command=gmx_run_cmd, task_work_path=ii, outlog='gmx_mdrun.log', errlog='gmx_mdrun.log') for ii in all_task_basedir ]
    gmx_run_submission = Submission(work_base=work_path, resources=resources, batch=machine, task_list=gmx_run_task)
    gmx_run_submission.run_submission()
