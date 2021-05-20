
import os, json, glob, shutil
import numpy as np
from rid.lib.utils import make_iter_name, make_walker_name, cmd_append_log, set_resource, set_machine, log_task

from rid.lib.cal_cv_dim import cal_cv_dim
from rid.lib.cmpf import cmpf
from rid.lib import LIB_PATH

from dpdispatcher.submission import Submission, Job, Task, Resources

enhc_out_conf="confs/"
enhc_out_angle="angle.rad.out"
enhc_out_plm="plm.out"
enhc_name="00.enhcMD"

res_name="01.resMD"
res_plm="plumed.res.dat"


def post_res (iter_index,
              json_file,
              cv_file,
              base_dir="./") :
    json_file = os.path.abspath(json_file)
    cv_file = os.path.abspath(cv_file)
    base_dir = os.path.abspath(base_dir) + "/"
    iter_name = make_iter_name (iter_index)
    res_path = base_dir + iter_name + "/" + res_name + "/"  
    cwd = os.getcwd()

    fp = open(json_file, 'r')
    jdata = json.load (fp)
    fp.close()

    os.chdir(res_path)
    all_task = glob.glob("/[0-9]*[0-9]")
    all_task = list(filter(lambda x:os.path.isdir(x),  glob.glob("[0-9]*[0-9]")))
    if len(all_task) == 0 :
        np.savetxt (res_path + 'data.raw', [], fmt = "%.6e")        
        os.chdir(cwd)
        return
    all_task.sort()
    all_task_reldir = [os.path.relpath(ii, res_path) for ii in all_task]
    
    centers = []
    force = []
    ndim = 0

    _conf_file = os.path.abspath(all_task[0] + "/conf.gro")
    cv_dim_list = cal_cv_dim (_conf_file, cv_file)
    cv_dih_dim = cv_dim_list[0]

    cmpf_cmd = "python3 {}/cmpf.py".format(LIB_PATH)
    cmpf_cmd += " -c %d" % cv_dih_dim
    cmpf_log = "cmpf.log"
    
    group_size = int((len(all_task)+1) // 8)
    print("rid.post_res.post_res:cmpf_cmd:", cmpf_cmd)
    print("rid.post_res.post_res:cmpf_group_size:", group_size)

    target = "cmpf"
    cmpf_resources =  Resources(
        number_node=jdata['{}_number_node'.format(target)], 
        cpu_per_node=jdata['{}_cpu_per_node'.format(target)], 
        gpu_per_node=jdata['{}_gpu_per_node'.format(target)], 
        queue_name=jdata['queue_name'], 
        group_size=group_size, 
        if_cuda_multi_devices=jdata['if_cuda_multi_devices']) 
    machine = set_machine(json_file)

    cmpf_task = [ Task(command=cmpf_cmd, task_work_path="{}".format(ii), outlog=cmpf_log, errlog=cmpf_log) for ii in all_task_reldir ]
    cmpf_submission = Submission(work_base=res_path, resources=cmpf_resources, batch=machine, task_list=cmpf_task)
    cmpf_submission.run_submission()
    print('cmpf done')

    abs_res_path = os.getcwd()
    for work_path in all_task:
        os.chdir (work_path)
        this_centers = np.loadtxt ('centers.out')
        centers = np.append (centers, this_centers)
        this_force = np.loadtxt ('force.out')
        force = np.append (force, this_force)        
        ndim = this_force.size
        assert (ndim == this_centers.size), "center size is diff to force size in " + work_path
        os.chdir(abs_res_path)

    os.chdir(cwd)
    centers = np.reshape (centers, [-1, ndim])
    force = np.reshape (force, [-1, ndim])
    data = np.concatenate ((centers, force), axis = 1)
    np.savetxt (res_path + 'data.raw', data, fmt = "%.6e")

    norm_force = np.linalg.norm (force, axis = 1)
    log_task ("min|f| = %e  max|f| = %e  avg|f| = %e" % 
              (np.min(norm_force), np.max(norm_force), np.average(norm_force)))
    print("min|f| = %e  max|f| = %e  avg|f| = %e" % (np.min(norm_force), np.max(norm_force), np.average(norm_force)))
    print('Saving cmpf finished!')
    print("Post process of restrained MD finished.")
    print(os.getcwd())


# def post_res_function (iter_index,
#               json_file,
#               cv_file,
#               base_dir="./") :
#     json_file = os.path.abspath(json_file)
#     cv_file = os.path.abspath(cv_file)
#     base_dir = os.path.abspath(base_dir) + "/"
#     iter_name = make_iter_name (iter_index)
#     res_path = base_dir + iter_name + "/" + res_name + "/"  
#     cwd = os.getcwd()

#     fp = open(json_file, 'r')
#     jdata = json.load (fp)
#     fp.close()

#     os.chdir(res_path)
#     all_task = glob.glob("/[0-9]*[0-9]")
#     all_task = list(filter(lambda x:os.path.isdir(x),  glob.glob("[0-9]*[0-9]")))
#     if len(all_task) == 0 :
#         np.savetxt (res_path + 'data.raw', [], fmt = "%.6e")        
#         return
#     all_task.sort()
    
#     cmpf_cmd = "python3 cmpf.py"
#     cmpf_log = "cmpf.log"
#     cmpf_cmd = cmd_append_log (cmpf_cmd, cmpf_log)

#     centers = []
#     force = []
#     ndim = 0
#     _conf_file = os.path.abspath(all_task[0] + "/conf.gro")
#     cv_dim_list = cal_cv_dim (_conf_file, cv_file)
#     cv_dih_dim = cv_dim_list[0]

#     for work_path in all_task:
#         os.chdir(work_path)
#         cmpf(cv_dih_dim, plm_out="plm.res.out", kappa_file='kappa.out', center_file='centers.out', tail=0.90, out_put='force.out')
#         this_centers = np.loadtxt ('centers.out')
#         centers = np.append (centers, this_centers)
#         this_force = np.loadtxt ('force.out')
#         force = np.append (force, this_force)        
#         ndim = this_force.size
#         assert (ndim == this_centers.size), "center size is diff to force size in " + work_path
#         os.chdir(res_path)
#     print('cmpf done')

#     os.chdir(cwd)
#     centers = np.reshape (centers, [-1, ndim])
#     force = np.reshape (force, [-1, ndim])
#     data = np.concatenate ((centers, force), axis = 1)
#     np.savetxt (res_path + 'data.raw', data, fmt = "%.6e")

#     norm_force = np.linalg.norm (force, axis = 1)
#     log_task ("min|f| = %e  max|f| = %e  avg|f| = %e" % 
#               (np.min(norm_force), np.max(norm_force), np.average(norm_force)))
#     print('Saving cmpf finished!')
#     print("Post process of restrained MD finished.")