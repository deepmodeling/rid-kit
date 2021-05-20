
import os, json, glob
from dpdispatcher.lazy_local_context import LazyLocalContext
from dpdispatcher.submission import Submission, Job, Task

from rid.lib.utils import make_iter_name, make_walker_name, cmd_append_log, set_resource, set_machine
from rid.lib.cal_cv_dim import cal_cv_dim
# from rid.lib.nn.train import train
# from rid.lib.nn.freeze import freeze_graph
from rid.lib.nn import NN_PATH

enhc_name="00.enhcMD"
train_name="02.train"


def check_new_data(iter_index, train_path, base_path):
    # check if new data is empty
    new_data_file = os.path.join(train_path, 'data/data.new.raw')
    if os.stat(new_data_file).st_size == 0 :
        prev_iter_index = iter_index - 1
        prev_train_path = base_path + make_iter_name(prev_iter_index) + "/" + train_name + "/"
        prev_models = glob.glob(prev_train_path + "*.pb")
        for ii in prev_models :
            model_name = os.path.basename(ii)
            os.symlink(ii, os.path.join(train_path, model_name))
        return True
    else:
        return False


def run_train (iter_index, 
               json_file, 
               cv_file,
               base_dir="./") :
    json_file = os.path.abspath(json_file)
    cv_file = os.path.abspath(cv_file)
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    fp.close()
    cv_file = os.path.abspath(cv_file)
    numb_model = jdata["numb_model"]
    train_thread = jdata["train_thread"]
    res_iter = jdata["res_iter"]
    base_dir = os.path.abspath(base_dir) + "/"
    iter_name = make_iter_name (iter_index)
    train_path = base_dir +  iter_name + "/" + train_name + "/"  
    if check_new_data(iter_index, train_path, base_dir):
        return

    enhc_path = base_dir + iter_name + "/" + enhc_name + "/"  
    _conf_file = enhc_path + "000/conf.gro"
    cv_dim_list = cal_cv_dim(_conf_file, cv_file)

    cwd = os.getcwd()
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
    cmdl_args += " -c "
    for cv_dim in cv_dim_list:
        cmdl_args += "%d " % cv_dim
    cmdl_args += " -b " + str(batch_size)
    cmdl_args += " -e " + str(numb_epoches)
    cmdl_args += " -l " + str(starter_lr)
    cmdl_args += " --decay-steps " + str(decay_steps)
    cmdl_args += " --decay-rate " + str(decay_rate)

    train_cmd = "python3 {}/train.py -t {:d}".format(NN_PATH, train_thread)
    train_cmd += cmdl_args
    train_cmd = cmd_append_log (train_cmd, "train.log")
    freez_cmd = "python3 {}/freeze.py -o graph.pb".format(NN_PATH)
    freez_cmd = cmd_append_log (freez_cmd, "freeze.log")    
    task_dirs = [ ("%03d" % ii) for ii in range(numb_model) ]
    
    print('lib.modeling.run_train:train_cmd:', train_cmd)
    print('lib.modeling.run_train:freez_cmd:', freez_cmd)
    print('lib.modeling.run_train:train_path:', train_path)
    print('lib.modeling.run_train:task_dirs:', task_dirs)
    
    # lazy_local_context = LazyLocalContext(local_root='./', work_profile=None)
    # pbs = PBS(context=lazy_local_context)
    # slurm = Slurm(context=lazy_local_context)
    resources = set_resource(json_file, target="train")
    machine = set_machine(json_file)

    train_task = [ Task(command=train_cmd, task_work_path=ii, outlog='train.log', errlog='train.log') for ii in task_dirs ]
    train_submission = Submission(work_base=train_path, resources=resources, batch=machine, task_list=train_task)
    train_submission.run_submission()

    freez_task = [ Task(command=freez_cmd, task_work_path=ii, outlog='freeze.log', errlog='freeze.log') for ii in task_dirs ]
    freez_submission = Submission(work_base=train_path, resources=resources, batch=machine, task_list=freez_task)
    freez_submission.run_submission()
    
    os.chdir(train_path)
    for ii in range(numb_model) :
        os.symlink ("%03d/graph.pb" % ii, "graph.%03d.pb" % ii)
    os.chdir(cwd)

    print("Training finished!")


# def run_train_function (iter_index, 
#                json_file, 
#                cv_file,
#                base_dir="./") :
#     json_file = os.path.abspath(json_file)
#     cv_file = os.path.abspath(cv_file)
#     fp = open (json_file, 'r')
#     jdata = json.load (fp)
#     fp.close()
#     cv_file = os.path.abspath(cv_file)
#     numb_model = jdata["numb_model"]
#     train_thread = jdata["train_thread"]
#     res_iter = jdata["res_iter"]
#     base_dir = os.path.abspath(base_dir) + "/"
#     iter_name = make_iter_name (iter_index)
#     train_path = base_dir +  iter_name + "/" + train_name + "/"  

#     enhc_path = base_dir + iter_name + "/" + enhc_name + "/"  
#     _conf_file = enhc_path + "000/conf.gro"
#     cwd = os.getcwd()
    
#     if check_new_data(iter_index, train_path, base_dir):
#         return
    
#     neurons = jdata["neurons"]
#     batch_size = jdata["batch_size"]
#     if iter_index < res_iter :
#         numb_epoches = jdata["numb_epoches"]
#         starter_lr = jdata["starter_lr"]
#         decay_steps = jdata["decay_steps"]
#         decay_rate = jdata["decay_rate"]    
#         cmdl_args = ""
#         use_mix = False
#         restart = False
#         old_ratio = 7.0
#     else :
#         numb_epoches = jdata["res_numb_epoches"]
#         starter_lr = jdata["res_starter_lr"]
#         decay_steps = jdata["res_decay_steps"]
#         decay_rate = jdata["res_decay_rate"]            
#         old_ratio = jdata["res_olddata_ratio"]
#         use_mix = True
#         restart = True

#     task_dirs = [ ("%03d" % ii) for ii in range(numb_model) ]
#     cv_dim_list = cal_cv_dim(_conf_file, cv_file)
    
#     os.chdir (train_path)
#     for work_path in task_dirs:
#         os.chdir(work_path)
#         train(
#             cv_dim_list,
#             neurons=neurons,
#             numb_threads=train_thread,
#             resnet=True,
#             use_mix=use_mix,
#             restart=restart,
#             batch_size=batch_size,
#             epoches=int(numb_epoches),
#             lr=float(starter_lr),
#             decay_steps=decay_steps,
#             decay_rate=decay_rate,
#             old_ratio=old_ratio,
#             decay_steps_inner=0,
#             init_model=None
#         )
#         freeze_graph("./", output="graph.pb")
#         os.chdir(train_path)

#     for ii in range(numb_model) :
#         os.symlink ("%03d/graph.pb" % ii, "graph.%03d.pb" % ii)
#     os.chdir (cwd)
#     print("Training finished!")