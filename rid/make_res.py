import os, json, glob, shutil
import numpy as np
from rid.lib.utils import create_path, replace, make_iter_name, make_walker_name, checkfile, log_task
from rid.lib.gen.gen_mdp import make_grompp
from rid.lib.gen.gen_plumed import make_res_templ_plumed, conf_res_plumed
from rid.lib.gen.gen_shell import gen_res_shell
from rid.lib.std import make_std
from rid.lib.cal_cv_dim import cal_cv_dim
from rid.lib.cluster_cv import sel_from_cluster

iter_format = "%06d"
walker_format = "%03d"
task_format = "%02d"
shell_clustering = False

enhc_name="00.enhcMD"
enhc_out_angle="angle.rad.out"
enhc_out_conf="confs/"
res_name="01.resMD"
res_plm="plumed.res.dat"

def make_sel_list(nconf, sel_idx):
    sel_list=""
    for ii in range (nconf) : 
        if ii == 0 : sel_list = str(sel_idx[ii])
        else : sel_list += "," + str(sel_idx[ii])
    return sel_list

def make_conf(nconf, res_path, walker_idx, walker_path, sel_idx, jdata, mol_path, conf_start=0, conf_every=1):
    mol_path = os.path.abspath(mol_path) + "/"
    mol_files = ["topol.top"]
    nsteps = jdata["res_nsteps"]
    frame_freq = jdata["res_frame_freq"]
    dt = jdata["res_dt"]
    temperature = jdata["res_temperature"]
    for ii in range (conf_start, nconf, conf_every) :
        work_path = res_path + ((walker_format + ".%06d") % (walker_idx, sel_idx[ii])) + "/"
        os.makedirs(work_path)
        make_grompp(work_path + "grompp.mdp", "res", nsteps, frame_freq, temperature=temperature, dt=dt, define="")
        for jj in mol_files:
            checkfile(work_path + jj)
            shutil.copy(mol_path + jj, work_path)

        conf_file = walker_path + enhc_out_conf + ("conf%d.gro" % sel_idx[int(ii)])
        checkfile(work_path + "conf.gro")
        tmp_cwd = os.getcwd()
        os.chdir(work_path)
        os.symlink (os.path.relpath(conf_file), "conf.gro")
        os.chdir(tmp_cwd)


def make_res_plumed(nconf, jdata, res_path, walker_idx, sel_idx, res_angles, conf_file, cv_file, conf_start=0, conf_every=1):
    res_kappa = jdata['res_kappa']
    res_ang_stride = jdata['res_ang_stride']
    res_prt_file = jdata['res_prt_file']
    frame_freq = jdata['res_frame_freq']
    cwd = os.getcwd()
    cv_file = os.path.abspath(cv_file)
    for ii in range (conf_start, nconf, conf_every) :
        work_path = os.path.abspath( res_path + ((walker_format + ".%06d") % (walker_idx, sel_idx[ii])) ) + "/"
        os.chdir(work_path)
        templ_plm_path = work_path + "plumed.res.templ"
        make_res_templ_plumed(templ_plm_path, conf_file, cv_file, res_kappa, res_ang_stride, res_prt_file)
        dir_str = ((walker_format + ".%06d") % (walker_idx, sel_idx[ii]))
        arg_str = np.array2string(res_angles[ii], 
                                    formatter={'float_kind':lambda x: "%.6f" % x}).replace("[","").replace("]","").replace("\n"," ")
        gen_res_shell("./general_mkres.sh")
        os.system("sh general_mkres.sh {}".format(arg_str))
        conf_res_plumed(work_path + res_plm, frame_freq)
        log_task (dir_str + ": " + arg_str)
        os.chdir(cwd)


def make_threshold(walker_idx, walker_path, base_path, sel_angles, cluster_threshold, init_numb_cluster, cv_dih_dim):
    if walker_idx ==0:
        cls_sel = sel_from_cluster (sel_angles, cluster_threshold, cv_dih_dim)
        test_numb_cluster=len(set(cls_sel))
        print('test number of clusters', test_numb_cluster)
        for test_iter in range(500):
            print('test_iter', test_iter)
            if test_numb_cluster < init_numb_cluster[0]:
                cluster_threshold=cluster_threshold * 0.95
                cls_sel = sel_from_cluster(sel_angles, cluster_threshold, cv_dih_dim)
                test_numb_cluster=len(set(cls_sel))
                print('cluster threshold', cluster_threshold)
                print('test number clusters', test_numb_cluster)
            elif test_numb_cluster > init_numb_cluster[1]:
                cluster_threshold=cluster_threshold * 1.05
                cls_sel = sel_from_cluster (sel_angles, cluster_threshold, cv_dih_dim)
                test_numb_cluster=len(set(cls_sel))
                print('cluster threshold', cluster_threshold)
                print('test number clusters', test_numb_cluster)
            else:
                print('cluster threshold', cluster_threshold)
                np.savetxt (walker_path + 'cluster_threshold.dat', [cluster_threshold], fmt = '%f')
                np.savetxt (base_path + 'cluster_threshold.dat', [cluster_threshold], fmt = '%f')
                break
        return cls_sel, cluster_threshold
    else:
        cluster_threshold = np.loadtxt(base_path + "cluster_threshold.dat")
        return None, cluster_threshold


def config_cls(sel_idx, cls_sel, max_sel, walker_path, cluster_threshold, sel_angles):
    if len(sel_idx) > 1:
        np.savetxt (walker_path + 'num_of_cluster.dat', [len(set(cls_sel))], fmt = '%d')
        np.savetxt (walker_path + 'cluster_threshold.dat', [cluster_threshold], fmt = '%f')
        if len(cls_sel)>max_sel:
            cls_sel=cls_sel[-max_sel:]
        sel_idx = sel_idx[cls_sel]
        np.savetxt (walker_path + 'cls.sel.angle.0.out', sel_angles[cls_sel], fmt = '%.6f')
    elif len(sel_idx)==1:
        np.savetxt (walker_path + 'num_of_cluster.dat', [1], fmt = '%d')
    return sel_idx


def make_res (iter_index, 
              json_file,
              cv_file,
              mol_path,
              base_dir="./") :

    json_file = os.path.abspath(json_file)
    fp = open (json_file, 'r')
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
    init_numb_cluster=[init_numb_cluster_lower, init_numb_cluster_upper]

    base_dir = os.path.abspath(base_dir) + "/"
    iter_name = make_iter_name (iter_index)
    enhc_path = base_dir + iter_name + "/" + enhc_name + "/" 
    res_path = base_dir + iter_name + "/" + res_name + "/"
    create_path (res_path) 

    cwd = os.getcwd()
    _conf_file = enhc_path + make_walker_name(0) + "/" + "conf.gro"
    cv_dim_list = cal_cv_dim (_conf_file, cv_file)
    cv_dim = sum(cv_dim_list)
    cv_dih_dim = cv_dim_list[0]
    ret_list = [True for ii in range(numb_walkers)]

    ## check if we have graph in enhc
    for walker_idx in range(numb_walkers) :
        cls_sel = None
        walker_path = enhc_path + walker_format % walker_idx + "/"
        graph_files = glob.glob (walker_path + "/*.pb")
        if len (graph_files) != 0 :
            cluster_threshold = np.loadtxt(base_dir + "cluster_threshold.dat")
            os.chdir (walker_path)
            models = glob.glob("*.pb")
            std_message = make_std(cv_dim, dataset=enhc_out_angle, models=models, threshold=sel_threshold, output="sel.out", output_angle="sel.angle.out")
            os.system('echo "{}" > sel.log'.format(std_message))
            log_task ("select with threshold %f" % sel_threshold)
            os.chdir (cwd)

            sel_idx = []
            sel_angles = np.array([])
            with open (walker_path + "sel.out") as fp :
                for line in fp : 
                    sel_idx += [int(x) for x in line.split()]
            if len(sel_idx) != 0 :
                sel_angles = np.reshape(np.loadtxt(walker_path + 'sel.angle.out'), [-1,cv_dim])
            elif len(sel_idx) == 0 :
                np.savetxt (walker_path + 'num_of_cluster.dat', [0], fmt = '%d')
                np.savetxt (walker_path + 'cls.sel.out', [], fmt = '%d')
                continue
        else :
            cluster_threshold = jdata["cluster_threshold"]
            sel_idx = range (len(glob.glob (walker_path + enhc_out_conf + "conf*gro")))
            sel_angles = np.loadtxt (walker_path + enhc_out_angle)
            sel_angles = np.reshape (sel_angles, [-1, cv_dim])            
            np.savetxt(walker_path + 'sel.out', sel_idx, fmt = '%d')
            np.savetxt(walker_path + 'sel.angle.out', sel_angles, fmt = '%.6f')
            cls_sel, cluster_threshold = make_threshold(walker_idx, walker_path, base_dir, sel_angles, cluster_threshold, init_numb_cluster, cv_dih_dim)
        if cls_sel is None:
            cls_sel = sel_from_cluster (sel_angles, cluster_threshold, cv_dih_dim)
        

        conf_start = 0
        conf_every = 1

        sel_idx = np.array (sel_idx, dtype = np.int)
        assert (len(sel_idx) == sel_angles.shape[0])
        sel_idx = config_cls(sel_idx, cls_sel, max_sel, walker_path, cluster_threshold, sel_angles)

        res_angles = np.loadtxt (walker_path + enhc_out_angle)
        res_angles = np.reshape (res_angles, [-1, cv_dim])
        res_angles = res_angles[sel_idx]
        np.savetxt (walker_path + 'cls.sel.out', sel_idx, fmt = '%d')
        np.savetxt (walker_path + 'cls.sel.angle.out', res_angles, fmt = '%.6f')
        res_confs = []
        for ii in sel_idx : 
            res_confs.append (walker_path + enhc_out_conf + ("conf%d.gro" % ii))    

        assert (len(res_confs) == res_angles.shape[0]), "number of enhc out conf does not match out angle"
        assert (len(sel_idx) == res_angles.shape[0]), "number of enhc out conf does not match number sel"
        nconf = len(res_confs)
        if nconf == 0 : 
            ret_list[walker_idx] = False
            continue

        sel_list = make_sel_list(nconf, sel_idx)
        log_task ("selected %d confs, indexes: %s" % (nconf, sel_list))
        make_conf(nconf, res_path, walker_idx, walker_path, sel_idx, jdata, mol_path, conf_start=0, conf_every=1)
        make_res_plumed(nconf, jdata, res_path, walker_idx, sel_idx, res_angles, _conf_file, cv_file, conf_start=0, conf_every=1)
    print("Restrained MD has prepared.")
