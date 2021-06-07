import os
import json
import glob
import shutil
from rid.lib.utils import create_path
from rid.lib.utils import log_task
from rid.lib.utils import make_iter_name, make_walker_name, cmd_append_log
from rid.lib.cal_cv_dim import cal_cv_dim
from rid.lib.nn import NN_PATH

from rid.lib.machine import set_resource, set_machine
from dpdispatcher.submission import Submission, Job, Task

enhc_name = "00.enhcMD"
res_name = "01.resMD"
train_name = "02.train"


def collect_data(iter_index, base_dir):
    iter_name = make_iter_name(iter_index)
    train_path = base_dir + iter_name + "/" + train_name + "/"
    data_path = train_path + "data/"
    data_file = train_path + "data/data.raw"
    data_old_file = train_path + "data/data.old.raw"
    data_new_file = train_path + "data/data.new.raw"
    cwd = os.getcwd() + "/"
    # collect data
    log_task("collect data upto %d" % (iter_index))
    if iter_index == 0:
        ii = 0
        this_raw = base_dir + make_iter_name(ii) + "/" + res_name + "/data.raw"
        os.chdir(data_path)
        os.symlink(os.path.relpath(this_raw), os.path.basename(data_new_file))
        os.symlink(os.path.basename(data_new_file),
                   os.path.basename(data_file))
        os.chdir(cwd)
        open(data_old_file, "w").close()
    else:
        prev_iter_index = iter_index - 1
        prev_data_file = base_dir + \
            make_iter_name(prev_iter_index) + "/" + \
            train_name + "/data/data.raw"
        this_raw = base_dir + \
            make_iter_name(iter_index) + "/" + res_name + "/data.raw"
        os.chdir(data_path)
        os.symlink(os.path.relpath(prev_data_file),
                   os.path.basename(data_old_file))
        os.symlink(os.path.relpath(this_raw), os.path.basename(data_new_file))
        os.chdir(cwd)
        with open(data_file, "wb") as fo:
            with open(data_old_file, "rb") as f0, open(data_new_file, "rb") as f1:
                shutil.copyfileobj(f0, fo)
                shutil.copyfileobj(f1, fo)


def make_train(iter_index,
               json_file,
               base_dir="./"):
    json_file = os.path.abspath(json_file)
    fp = open(json_file, 'r')
    jdata = json.load(fp)
    fp.close()
    # template_dir = jdata["template_dir"]
    numb_model = jdata["numb_model"]
    res_iter = jdata["res_iter"]

    # abs path
    base_dir = os.path.abspath(base_dir) + "/"
    iter_name = make_iter_name(iter_index)
    train_path = base_dir + iter_name + "/" + train_name + "/"
    data_path = train_path + "data/"
    cwd = os.getcwd() + "/"
    create_path(train_path)
    os.makedirs(data_path)
    collect_data(iter_index, base_dir)

    # create train dirs
    log_task("create train dirs")
    for ii in range(numb_model):
        work_path = train_path + ("%03d/" % ii)
        old_model_path = work_path + "old_model/"
        create_path(work_path)
        os.chdir(work_path)
        os.symlink("../data", "./data")
        os.chdir(cwd)
        if iter_index >= 1:
            prev_iter_index = iter_index - 1
            prev_iter_name = make_iter_name(prev_iter_index)
            prev_train_path = base_dir + prev_iter_name + "/" + train_name + "/"
            prev_work_path = prev_train_path + ("%03d/" % ii)
            prev_model_files = glob.glob(prev_work_path + "model.ckpt.*")
            prev_model_files = prev_model_files + \
                [prev_work_path + "checkpoint"]
            create_path(old_model_path)
            os.chdir(old_model_path)
            # why to copy twice.
            for ii in prev_model_files:
                os.symlink(os.path.relpath(ii), os.path.basename(ii))
            os.chdir(cwd)
            for jj in prev_model_files:
                shutil.copy(jj, work_path)
    print("Training files have prepared.")


def check_new_data(iter_index, train_path, base_path):
    # check if new data is empty
    new_data_file = os.path.join(train_path, 'data/data.new.raw')
    if os.stat(new_data_file).st_size == 0:
        prev_iter_index = iter_index - 1
        prev_train_path = base_path + \
            make_iter_name(prev_iter_index) + "/" + train_name + "/"
        prev_models = glob.glob(prev_train_path + "*.pb")
        for ii in prev_models:
            model_name = os.path.basename(ii)
            os.symlink(ii, os.path.join(train_path, model_name))
        return True
    else:
        return False


def run_train(iter_index,
              json_file,
              machine_json,
              cv_file,
              base_dir="./"):
    json_file = os.path.abspath(json_file)
    cv_file = os.path.abspath(cv_file)
    fp = open(json_file, 'r')
    jdata = json.load(fp)
    fp.close()
    cv_file = os.path.abspath(cv_file)
    numb_model = jdata["numb_model"]
    train_thread = jdata["train_thread"]
    res_iter = jdata["res_iter"]
    base_dir = os.path.abspath(base_dir) + "/"
    iter_name = make_iter_name(iter_index)
    train_path = base_dir + iter_name + "/" + train_name + "/"
    if check_new_data(iter_index, train_path, base_dir):
        return

    enhc_path = base_dir + iter_name + "/" + enhc_name + "/"
    _conf_file = enhc_path + "000/conf.gro"
    cv_dim_list = cal_cv_dim(_conf_file, cv_file)

    cwd = os.getcwd()
    neurons = jdata["neurons"]
    batch_size = jdata["batch_size"]
    if iter_index < res_iter:
        numb_epoches = jdata["numb_epoches"]
        starter_lr = jdata["starter_lr"]
        decay_steps = jdata["decay_steps"]
        decay_rate = jdata["decay_rate"]
        cmdl_args = ""
    else:
        numb_epoches = jdata["res_numb_epoches"]
        starter_lr = jdata["res_starter_lr"]
        decay_steps = jdata["res_decay_steps"]
        decay_rate = jdata["res_decay_rate"]
        old_ratio = jdata["res_olddata_ratio"]
        cmdl_args = " --restart --use-mix --old-ratio %f " % old_ratio

    if jdata["resnet"]:
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
    train_cmd = cmd_append_log(train_cmd, "train.log")
    freez_cmd = "python3 {}/freeze.py -o graph.pb".format(NN_PATH)
    freez_cmd = cmd_append_log(freez_cmd, "freeze.log")
    task_dirs = [("%03d" % ii) for ii in range(numb_model)]

    print('lib.modeling.run_train:train_cmd:', train_cmd)
    print('lib.modeling.run_train:freez_cmd:', freez_cmd)
    print('lib.modeling.run_train:train_path:', train_path)
    print('lib.modeling.run_train:task_dirs:', task_dirs)

    resources = set_resource(machine_json, target="train")
    machine = set_machine(machine_json, target="train")

    train_task = [Task(command=train_cmd, task_work_path=ii,
                       outlog='train.log', errlog='train.log') for ii in task_dirs]
    train_submission = Submission(
        work_base=train_path, machine=machine, resources=resources, task_list=train_task)
    train_submission.run_submission()

    freez_task = [Task(command=freez_cmd, task_work_path=ii,
                       outlog='freeze.log', errlog='freeze.log') for ii in task_dirs]
    freez_submission = Submission(
        work_base=train_path, machine=machine, resources=resources, task_list=freez_task)
    freez_submission.run_submission()

    os.chdir(train_path)
    for ii in range(numb_model):
        os.symlink("%03d/graph.pb" % ii, "graph.%03d.pb" % ii)
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
