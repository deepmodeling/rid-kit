import os
import json
import glob
import shutil
import numpy as np
from rid.lib.utils import create_path, replace, make_iter_name, make_walker_name, checkfile, log_task, cmd_append_log
from rid.lib.gen.gen_mdp import make_grompp
from rid.lib.gen.gen_plumed import make_res_templ_plumed, conf_res_plumed
from rid.lib.gen.gen_shell import gen_res_shell
from rid.lib.std import make_std
from rid.lib.cal_cv_dim import cal_cv_dim
from rid.lib.cluster_cv import sel_from_cluster

from rid.lib.cmpf import cmpf
from rid.lib import LIB_PATH

from rid.lib.machine import set_resource, set_machine
from dpdispatcher.submission import Submission, Job, Task

res_plm = "plumed.res.dat"
res_name = "01.resMD"

enhc_out_conf = "confs/"
enhc_out_angle = "angle.rad.out"
enhc_out_plm = "plm.out"
enhc_name = "00.enhcMD"

iter_format = "%06d"
walker_format = "%03d"
task_format = "%02d"
shell_clustering = False


def make_sel_list(nconf, sel_idx):
    sel_list = ""
    for ii in range(nconf):
        if ii == 0:
            sel_list = str(sel_idx[ii])
        else:
            sel_list += "," + str(sel_idx[ii])
    return sel_list


def make_conf(nconf, res_path, walker_idx, walker_path, sel_idx, jdata, mol_path, conf_start=0, conf_every=1):
    mol_path = os.path.abspath(mol_path) + "/"
    mol_files = ["topol.top"]
    nsteps = jdata["res_nsteps"]
    frame_freq = jdata["res_frame_freq"]
    dt = jdata["res_dt"]
    temperature = jdata["res_temperature"]
    for ii in range(conf_start, nconf, conf_every):
        work_path = res_path + ((walker_format + ".%06d") %
                                (walker_idx, sel_idx[ii])) + "/"
        os.makedirs(work_path)
        make_grompp(work_path + "grompp.mdp", "res", nsteps,
                    frame_freq, temperature=temperature, dt=dt, define="")
        for jj in mol_files:
            checkfile(work_path + jj)
            shutil.copy(mol_path + jj, work_path)

        conf_file = walker_path + enhc_out_conf + \
            ("conf%d.gro" % sel_idx[int(ii)])
        checkfile(work_path + "conf.gro")
        tmp_cwd = os.getcwd()
        os.chdir(work_path)
        os.symlink(os.path.relpath(conf_file), "conf.gro")
        os.chdir(tmp_cwd)


def make_res_plumed(nconf, jdata, res_path, walker_idx, sel_idx, res_angles, conf_file, cv_file, conf_start=0, conf_every=1):
    res_kappa = jdata['res_kappa']
    res_ang_stride = jdata['res_ang_stride']
    res_prt_file = jdata['res_prt_file']
    frame_freq = jdata['res_frame_freq']
    cwd = os.getcwd()
    cv_file = os.path.abspath(cv_file)
    for ii in range(conf_start, nconf, conf_every):
        work_path = os.path.abspath(
            res_path + ((walker_format + ".%06d") % (walker_idx, sel_idx[ii]))) + "/"
        os.chdir(work_path)
        templ_plm_path = work_path + "plumed.res.templ"
        make_res_templ_plumed(templ_plm_path, conf_file,
                              cv_file, res_kappa, res_ang_stride, res_prt_file)
        dir_str = ((walker_format + ".%06d") % (walker_idx, sel_idx[ii]))
        arg_str = np.array2string(res_angles[ii],
                                  formatter={'float_kind': lambda x: "%.6f" % x}).replace("[", "").replace("]", "").replace("\n", " ")
        gen_res_shell("./general_mkres.sh")
        os.system("sh general_mkres.sh {}".format(arg_str))
        conf_res_plumed(work_path + res_plm, frame_freq)
        log_task(dir_str + ": " + arg_str)
        os.chdir(cwd)


def make_threshold(walker_idx, walker_path, base_path, sel_angles, cluster_threshold, init_numb_cluster, cv_dih_dim, weight):
    if walker_idx == 0:
        cls_sel = sel_from_cluster(sel_angles, cluster_threshold, cv_dih_dim, weight)
        test_numb_cluster = len(set(cls_sel))
        print('test number of clusters', test_numb_cluster)
        for test_iter in range(500):
            print('test_iter', test_iter)
            if test_numb_cluster < init_numb_cluster[0]:
                cluster_threshold = cluster_threshold * 0.95
                cls_sel = sel_from_cluster(
                    sel_angles, cluster_threshold, cv_dih_dim, weight)
                test_numb_cluster = len(set(cls_sel))
                print('cluster threshold', cluster_threshold)
                print('test number clusters', test_numb_cluster)
            elif test_numb_cluster > init_numb_cluster[1]:
                cluster_threshold = cluster_threshold * 1.05
                cls_sel = sel_from_cluster(
                    sel_angles, cluster_threshold, cv_dih_dim, weight)
                test_numb_cluster = len(set(cls_sel))
                print('cluster threshold', cluster_threshold)
                print('test number clusters', test_numb_cluster)
            else:
                print('cluster threshold', cluster_threshold)
                np.savetxt(walker_path + 'cluster_threshold.dat',
                           [cluster_threshold], fmt='%f')
                np.savetxt(base_path + 'cluster_threshold.dat',
                           [cluster_threshold], fmt='%f')
                break
        return cls_sel, cluster_threshold
    else:
        cluster_threshold = np.loadtxt(base_path + "cluster_threshold.dat")
        return None, cluster_threshold


def config_cls(sel_idx, cls_sel, max_sel, walker_path, cluster_threshold, sel_angles):
    if len(sel_idx) > 1:
        np.savetxt(walker_path + 'num_of_cluster.dat',
                   [len(set(cls_sel))], fmt='%d')
        np.savetxt(walker_path + 'cluster_threshold.dat',
                   [cluster_threshold], fmt='%f')
        if len(cls_sel) > max_sel:
            cls_sel = cls_sel[-max_sel:]
        sel_idx = sel_idx[cls_sel]
        np.savetxt(walker_path + 'cls.sel.angle.0.out',
                   sel_angles[cls_sel], fmt='%.6f')
    elif len(sel_idx) == 1:
        np.savetxt(walker_path + 'num_of_cluster.dat', [1], fmt='%d')
    return sel_idx


def make_res(iter_index,
             json_file,
             cv_file,
             mol_path,
             base_dir="./"):

    json_file = os.path.abspath(json_file)
    fp = open(json_file, 'r')
    jdata = json.load(fp)
    fp.close()
    cv_file = os.path.abspath(cv_file)

    numb_walkers = jdata["numb_walkers"]
    bias_nsteps = jdata["bias_nsteps"]
    bias_frame_freq = jdata["bias_frame_freq"]
    nsteps = jdata["res_nsteps"]
    frame_freq = jdata["res_frame_freq"]
    sel_threshold = jdata["sel_threshold"]
    max_sel = jdata["max_sel"]
    cluster_threshold = jdata["cluster_threshold"]
    init_numb_cluster_upper = int(jdata["init_numb_cluster_upper"])
    init_numb_cluster_lower = int(jdata["init_numb_cluster_lower"])
    init_numb_cluster = [init_numb_cluster_lower, init_numb_cluster_upper]

    base_dir = os.path.abspath(base_dir) + "/"
    iter_name = make_iter_name(iter_index)
    enhc_path = base_dir + iter_name + "/" + enhc_name + "/"
    res_path = base_dir + iter_name + "/" + res_name + "/"
    create_path(res_path)

    cwd = os.getcwd()
    _conf_file = enhc_path + make_walker_name(0) + "/" + "conf.gro"
    cv_dim_list = cal_cv_dim(_conf_file, cv_file)
    cv_dim = sum(cv_dim_list)
    cv_dih_dim = cv_dim_list[0]
    ret_list = [True for ii in range(numb_walkers)]

    weight = jdata["cv_weight_for_cluster"]
    if type(weight) == list:
        assert len(weight) == cv_dim, "Number of values in the weight list is not equal to the number of CVs."
    elif type(weight) == float or type(weight) == int:
        assert weight != 0
    else:
        raise TypeError("Invalid type of weight of CVs for clustering. Please use int or list instead.")

    # check if we have graph in enhc
    for walker_idx in range(numb_walkers):
        cls_sel = None
        walker_path = enhc_path + walker_format % walker_idx + "/"
        graph_files = glob.glob(walker_path + "*.pb")
        if len(graph_files) != 0:
            cluster_threshold = np.loadtxt(base_dir + "cluster_threshold.dat")
            os.chdir(walker_path)
            models = glob.glob("*.pb")
            std_message = make_std(cv_dim, dataset=enhc_out_angle, models=models,
                                   threshold=sel_threshold, output="sel.out", output_angle="sel.angle.out")
            os.system('echo "{}" > sel.log'.format(std_message))
            log_task("select with threshold %f" % sel_threshold)
            os.chdir(cwd)

            sel_idx = []
            sel_angles = np.array([])
            with open(walker_path + "sel.out") as fp:
                for line in fp:
                    sel_idx += [int(x) for x in line.split()]
            if len(sel_idx) != 0:
                sel_angles = np.reshape(np.loadtxt(
                    walker_path + 'sel.angle.out'), [-1, cv_dim])
            elif len(sel_idx) == 0:
                np.savetxt(walker_path + 'num_of_cluster.dat', [0], fmt='%d')
                np.savetxt(walker_path + 'cls.sel.out', [], fmt='%d')
                continue
        else:
            cluster_threshold = jdata["cluster_threshold"]
            sel_idx = range(
                len(glob.glob(walker_path + enhc_out_conf + "conf*gro")))
            sel_angles = np.loadtxt(walker_path + enhc_out_angle)
            sel_angles = np.reshape(sel_angles, [-1, cv_dim])
            np.savetxt(walker_path + 'sel.out', sel_idx, fmt='%d')
            np.savetxt(walker_path + 'sel.angle.out', sel_angles, fmt='%.6f')
            cls_sel, cluster_threshold = make_threshold(
                walker_idx, walker_path, base_dir, sel_angles, cluster_threshold, init_numb_cluster, cv_dih_dim, weight)
        if cls_sel is None:
            print(sel_angles, cluster_threshold, cv_dih_dim)
            cls_sel = sel_from_cluster(
                sel_angles, cluster_threshold, cv_dih_dim, weight)

        conf_start = 0
        conf_every = 1

        sel_idx = np.array(sel_idx, dtype=np.int)
        assert (len(sel_idx) == sel_angles.shape[0]), "{} selected indexes don't match {} selected angles.".format(len(sel_idx), sel_angles.shape[0])
        sel_idx = config_cls(sel_idx, cls_sel, max_sel,
                             walker_path, cluster_threshold, sel_angles)

        res_angles = np.loadtxt(walker_path + enhc_out_angle)
        res_angles = np.reshape(res_angles, [-1, cv_dim])
        res_angles = res_angles[sel_idx]
        np.savetxt(walker_path + 'cls.sel.out', sel_idx, fmt='%d')
        np.savetxt(walker_path + 'cls.sel.angle.out', res_angles, fmt='%.6f')
        res_confs = []
        for ii in sel_idx:
            res_confs.append(walker_path + enhc_out_conf + ("conf%d.gro" % ii))

        assert (len(res_confs) ==
                res_angles.shape[0]), "number of enhc out conf does not match out angle"
        assert (len(sel_idx) ==
                res_angles.shape[0]), "number of enhc out conf does not match number sel"
        nconf = len(res_confs)
        if nconf == 0:
            ret_list[walker_idx] = False
            continue

        sel_list = make_sel_list(nconf, sel_idx)
        log_task("selected %d confs, indexes: %s" % (nconf, sel_list))
        make_conf(nconf, res_path, walker_idx, walker_path, sel_idx,
                  jdata, mol_path, conf_start=0, conf_every=1)
        make_res_plumed(nconf, jdata, res_path, walker_idx, sel_idx,
                        res_angles, _conf_file, cv_file, conf_start=0, conf_every=1)
    print("Restrained MD has been prepared.")


def run_res(iter_index,
            json_file,
            machine_json,
            base_dir="./"):
    json_file = os.path.abspath(json_file)
    fp = open(json_file, 'r')
    jdata = json.load(fp)
    fp.close()
    gmx_prep = jdata["gmx_prep"]
    gmx_run = jdata["gmx_run"]
    res_thread = jdata["res_thread"]
    gmx_run = gmx_run + (" -nt %d" % res_thread)
    gmx_run = gmx_run + " -plumed " + res_plm
    # gmx_cont_run = gmx_run + " -cpi state.cpt"
    gmx_prep_log = "gmx_grompp.log"
    gmx_run_log = "gmx_mdrun.log"
    gmx_prep_cmd = cmd_append_log(gmx_prep, gmx_prep_log)
    gmx_run_cmd = cmd_append_log(gmx_run, gmx_run_log)
    # gmx_cont_run_cmd = cmd_append_log (gmx_cont_run, gmx_run_log)

    base_dir = os.path.abspath(base_dir) + "/"
    iter_name = make_iter_name(iter_index)
    res_path = base_dir + iter_name + "/" + res_name + "/"

    if not os.path.isdir(res_path):
        raise RuntimeError(
            "do not see any restrained simulation (%s)." % res_path)

    all_task = list(filter(lambda x: os.path.isdir(
        x),  glob.glob(res_path + "/[0-9]*[0-9]")))
    print('run_res:all_task_propose:', all_task)
    print('run_res:gmx_prep_cmd:', gmx_prep_cmd)
    print('run_res:gmx_run_cmd:', gmx_run_cmd)
    # print('run_res:gmx_cont_run_cmd:', gmx_cont_run_cmd)

    if len(all_task) == 0:
        return None
    all_task.sort()
    all_task_basedir = [os.path.relpath(ii, res_path) for ii in all_task]

    res_resources = set_resource(machine_json, target="resMD")
    machine = set_machine(machine_json, target="resMD")

    gmx_prep_task = [Task(command=gmx_prep_cmd, task_work_path=ii,
                          outlog='gmx_grompp.log', errlog='gmx_grompp.log') for ii in all_task_basedir]
    gmx_prep_submission = Submission(
        work_base=res_path, machine=machine, resources=res_resources, task_list=gmx_prep_task)
    gmx_prep_submission.run_submission()

    gmx_run_task = [Task(command=gmx_run_cmd, task_work_path=ii,
                         outlog='gmx_mdrun.log', errlog='gmx_mdrun.log') for ii in all_task_basedir]
    gmx_run_submission = Submission(
        work_base=res_path, machine=machine, resources=res_resources, task_list=gmx_run_task)
    gmx_run_submission.run_submission()


def post_res(iter_index,
             json_file,
             machine_json,
             cv_file,
             base_dir="./"):
    json_file = os.path.abspath(json_file)
    machine_json = os.path.abspath(machine_json)
    cv_file = os.path.abspath(cv_file)
    base_dir = os.path.abspath(base_dir) + "/"
    iter_name = make_iter_name(iter_index)
    res_path = base_dir + iter_name + "/" + res_name + "/"
    cwd = os.getcwd()

    fp = open(json_file, 'r')
    jdata = json.load(fp)
    fp.close()

    os.chdir(res_path)
    all_task = glob.glob("/[0-9]*[0-9]")
    all_task = list(filter(lambda x: os.path.isdir(x),
                           glob.glob("[0-9]*[0-9]")))
    if len(all_task) == 0:
        np.savetxt(res_path + 'data.raw', [], fmt="%.6e")
        os.chdir(cwd)
        return
    all_task.sort()
    all_task_reldir = [os.path.relpath(ii, res_path) for ii in all_task]

    centers = []
    force = []
    ndim = 0

    _conf_file = os.path.abspath(all_task[0] + "/conf.gro")
    cv_dim_list = cal_cv_dim(_conf_file, cv_file)
    cv_dih_dim = cv_dim_list[0]

    cmpf_cmd = "python3 {}/cmpf.py".format(LIB_PATH)
    cmpf_cmd += " -c %d" % cv_dih_dim
    cmpf_log = "cmpf.log"

    print("rid.post_res.post_res:cmpf_cmd:", cmpf_cmd)

    cmpf_resources = set_resource(machine_json, target="cmpf")
    machine = set_machine(machine_json, target="cmpf")

    cmpf_task = [Task(command=cmpf_cmd, task_work_path="{}".format(
        ii), outlog=cmpf_log, errlog=cmpf_log) for ii in all_task_reldir]
    cmpf_submission = Submission(
        work_base=res_path, machine=machine, resources=cmpf_resources, task_list=cmpf_task)
    cmpf_submission.run_submission()
    print('cmpf done')

    abs_res_path = os.getcwd()
    for work_path in all_task:
        os.chdir(work_path)
        this_centers = np.loadtxt('centers.out')
        centers = np.append(centers, this_centers)
        this_force = np.loadtxt('force.out')
        force = np.append(force, this_force)
        ndim = this_force.size
        assert (
            ndim == this_centers.size), "center size is diff to force size in " + work_path
        os.chdir(abs_res_path)

    os.chdir(cwd)
    centers = np.reshape(centers, [-1, ndim])
    force = np.reshape(force, [-1, ndim])
    data = np.concatenate((centers, force), axis=1)
    np.savetxt(res_path + 'data.raw', data, fmt="%.6e")

    norm_force = np.linalg.norm(force, axis=1)
    log_task("min|f| = %e  max|f| = %e  avg|f| = %e" %
             (np.min(norm_force), np.max(norm_force), np.average(norm_force)))
    print("min|f| = %e  max|f| = %e  avg|f| = %e" %
          (np.min(norm_force), np.max(norm_force), np.average(norm_force)))
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
