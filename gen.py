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
from tools.general_plumed import cal_cv_dim
from tools.general_plumed import general_plumed

mol_dir_fmt = "mol/alan/%s/%02d"
bf_file_copy = ["conf.gro", "grompp.mdp", "topol.top"]
res_file_copy = ["clean.sh", "cmpf.sh", "cmpf.py", "general_mkres.sh", "tools", "run.res.sh", "equi.mdp"]
afed_file_copy = ["sel.angle.py"]
rid_run = "rid.py"
rid_param = "rid.json"
rid_file_copy = [rid_run, rid_param, "template", "lib"]
rit_run = "rit.py"
rit_param = "rit.json"
rit_file_copy = [rit_run, rit_param, "template", "lib"]
strm_file_copy = ["strm.py", "init.py", "ener.py", "param.json", "librun.py", "libstring.py", "MachineSlurm.py", "MachineLocal.py"]

def replace (file_name, pattern, subst) :
    file_handel = open (file_name, 'r')
    file_string = file_handel.read ()
    file_handel.close ()
    file_string = ( re.sub (pattern, subst, file_string) )
    file_handel = open (file_name, 'w')
    file_handel.write (file_string)
    file_handel.close ()

def make_grompp_bias (gro_file, nsteps, frame_freq) :
    replace (gro_file, "nsteps.*=.*", "nsteps = %d" % nsteps)
    replace (gro_file, "nstxout.*=.*", "nstxout = %d" % frame_freq)
    replace (gro_file, "nstvout.*=.*", "nstvout = %d" % frame_freq)
    replace (gro_file, "nstfout.*=.*", "nstfout = %d" % frame_freq)
    replace (gro_file, "nstxtcout.*=.*", "nstxtcout = %d" % frame_freq)
    replace (gro_file, "nstenergy.*=.*", "nstenergy = %d" % frame_freq)    

def make_grompp_res (gro_file, nsteps, frame_freq) :
    replace (gro_file, "nsteps.*=.*", "nsteps = %d" % nsteps)
    replace (gro_file, "nstxout.*=.*", "nstxout = %d" % 0)
    replace (gro_file, "nstvout.*=.*", "nstvout = %d" % 0)
    replace (gro_file, "nstfout.*=.*", "nstfout = %d" % 0)
    replace (gro_file, "nstxtcout.*=.*", "nstxtcout = %d" % frame_freq)
    replace (gro_file, "nstenergy.*=.*", "nstenergy = %d" % frame_freq)

def create_path (path) :
    if os.path.isdir(path) : 
        dirname = os.path.dirname(path)
        counter = 0
        while True :
            bk_dirname = dirname + ".bk%03d" % counter
            if not os.path.isdir(bk_dirname) : 
                shutil.move (dirname, bk_dirname) 
                break
            counter += 1
    os.makedirs (path)

def gen_bf (out_dir, 
            gen_file, 
            cv_file, 
            mol_dir) :    
    fp = open (gen_file, 'r')
    jdata = json.load (fp)    
    bf_traj_stride = jdata['bf_traj_stride']

    mol_dir += "/"
    out_dir += "/"
    
    create_path (out_dir)

    for ii in bf_file_copy :        
        shutil.copy (mol_dir + ii, out_dir)
    make_grompp_bias (out_dir + "grompp.mdp", 50000, bf_traj_stride)

    confs = glob.glob (mol_dir + "conf*gro")
    for cc in confs:
        shutil.copy (cc, out_dir)    

    ret = general_plumed ("bf", confs[0], cv_file, pstride = bf_traj_stride)
    with open(out_dir + "plumed.dat", "w") as fp:
        fp.write (ret)

    # angle_idxes_arg=""
    # angle_idxes = np.reshape(angle_idxes, [-1])
    # for ii in angle_idxes :
    #     angle_idxes_arg += " " + str(ii)
    # sp.check_call (tools_dir + "gen.anglendx.py " + 
    #                " -i " + angle_idxes_arg +
    #                " -a " + str(nalpha) + 
    #                " --alpha-fmt " + alpha_idx_fmt + 
    #                " --angle-fmt " + angle_idx_fmt + 
    #                " > " + out_dir + "angle.ndx",
    #                shell = True
    # )    
    # shutil.copy (tools_dir + "angle.sh", out_dir)    
    # replace (out_dir + "angle.sh", "nalpha=.*", "nalpha=%d" % nalpha)
    # replace (out_dir + "angle.sh", "nangle=.*", "nangle=%d" % nangle)

def gen_res (out_dir,
             gen_file,
             cv_file,
             mol_dir,
             res_dir) :
    fp = open (gen_file, 'r')
    jdata = json.load (fp)    
    res_kappa = jdata['res_kappa']
    res_traj_stride = jdata['res_traj_stride']
    res_ang_stride = jdata['res_ang_stride']
    res_prt_file = jdata['res_prt_file']

    mol_dir += "/"
    res_dir += "/"
    out_dir += "/"
    conf_file = mol_dir + "conf.gro"

    create_path (out_dir)
    
    for ii in bf_file_copy :        
        shutil.copy (mol_dir + ii, out_dir)
    make_grompp_res (out_dir + "grompp.mdp", 50000, res_traj_stride)

    for ii in res_file_copy :        
        if os.path.isdir (res_dir + ii) :
            shutil.copytree (res_dir + ii, out_dir + ii)
        elif os.path.isfile (res_dir + ii) :
            shutil.copy (res_dir + ii, out_dir)

    ret = general_plumed ("res", conf_file, cv_file, 
                          kappa = res_kappa,
                          pstride = res_ang_stride, 
                          pfile = res_prt_file)
    with open(out_dir + "plumed.res.templ", "w") as fp:
        fp.write (ret)

    conf_file = mol_dir + "conf.gro"
    cv_dim_list = cal_cv_dim (conf_file, cv_file)
    replace (out_dir + "/cmpf.py", "cv_dih_dim = .*", "cv_dih_dim = %d" % cv_dim_list[0])
    
    # replace (out_dir + "mkres.sh", "alpha_idx_fmt=.*", "alpha_idx_fmt=%s" % alpha_idx_fmt)
    # replace (out_dir + "mkres.sh", "angle_idx_fmt=.*", "angle_idx_fmt=%s" % angle_idx_fmt)
    # replace (out_dir + "mkres.sh", "nalpha=.*", "nalpha=%d" % nalpha)
    # replace (out_dir + "mkres.sh", "nangle=.*", "nangle=%d" % nangle)

    # angle_idxes_arg=""
    # angle_idxes = np.reshape(angle_idxes, [-1])
    # for ii in angle_idxes :
    #     angle_idxes_arg += " " + str(ii)    
    # sp.check_call (tools_dir + "gen.anglendx.py " + 
    #                " -i " + angle_idxes_arg +
    #                " -a " + str(nalpha) + 
    #                " --alpha-fmt " + alpha_idx_fmt + 
    #                " --angle-fmt " + angle_idx_fmt + 
    #                " > " + out_dir + "angle.ndx",
    #                shell = True
    # )
    # sp.check_call (tools_dir + "gen.plumed.py " +
    #                " res " +
    #                " -i " + angle_idxes_arg + 
    #                " -a " + str(nalpha) + 
    #                " --alpha-fmt " + alpha_idx_fmt + 
    #                " --angle-fmt " + angle_idx_fmt + 
    #                " -k " + str(res_kappa) + 
    #                " -s " + str(res_ang_stride) + 
    #                " -f " + str(res_prt_file) + 
    #                " > " + out_dir + "plumed.res.templ",
    #                shell = True
    # )
    # shutil.copy (tools_dir + "angle.sh", out_dir)    
    # replace (out_dir + "angle.sh", "nalpha=.*", "nalpha=%d" % nalpha)
    # replace (out_dir + "angle.sh", "nangle=.*", "nangle=%d" % nangle)    


def gen_afed (out_dir,
              gen_file,
              cv_file,
              mol_dir,
              afed_dir) :
    fp = open (gen_file, 'r')
    jdata = json.load (fp)    
    afed_temp = jdata['afed_temp']
    afed_kappa = jdata['afed_kappa']
    afed_gamma = jdata['afed_gamma']
    afed_tau = jdata['afed_tau']
    afed_traj_stride = jdata['afed_traj_stride']
    afed_prt_file = jdata['afed_prt_file']

    mol_dir += "/"
    afed_dir += "/"
    out_dir += "/"
    conf_file = mol_dir + "conf.gro"

    create_path (out_dir)
    
    for ii in bf_file_copy :        
        shutil.copy (mol_dir + ii, out_dir)
    make_grompp_bias (out_dir + "grompp.mdp", 50000, afed_traj_stride)

    for ii in afed_file_copy :        
        if os.path.isdir (afed_dir + ii) :
            shutil.copytree (afed_dir + ii, out_dir + ii)
        elif os.path.isfile (afed_dir + ii) :
            shutil.copy (afed_dir + ii, out_dir)

    ret = general_plumed("afed", conf_file, cv_file, 
                         kappa = afed_kappa,
                         temp = afed_temp,
                         tau = afed_tau,
                         gamma = afed_gamma,
                         pstride = afed_traj_stride,
                         pfile = afed_prt_file)
    with open(out_dir + "plumed.afed.dat", "w") as fp:
        fp.write (ret)
    
    # angle_idxes_arg=""
    # angle_idxes = np.reshape(angle_idxes, [-1])
    # for ii in angle_idxes :
    #     angle_idxes_arg += " " + str(ii)    
    # sp.check_call (tools_dir + "gen.anglendx.py " + 
    #                " -i " + angle_idxes_arg +
    #                " -a " + str(nalpha) + 
    #                " --alpha-fmt " + alpha_idx_fmt + 
    #                " --angle-fmt " + angle_idx_fmt + 
    #                " > " + out_dir + "angle.ndx",
    #                shell = True
    # )
    # sp.check_call (tools_dir + "gen.plumed.py " +
    #                " afed " +
    #                " -i " + angle_idxes_arg + 
    #                " -a " + str(nalpha) + 
    #                " --alpha-fmt " + alpha_idx_fmt + 
    #                " --angle-fmt " + angle_idx_fmt + 
    #                " -T " + str(afed_temp) + 
    #                " -k " + str(afed_kappa) + 
    #                " -g " + str(afed_gamma) + 
    #                " -t " + str(afed_tau) + 
    #                " -s " + str(afed_traj_stride) + 
    #                " -f " + str(afed_prt_file) + 
    #                " > " + out_dir + "plumed.afed.dat",
    #                shell = True
    # )
    # shutil.copy (tools_dir + "angle.sh", out_dir)
    # replace (out_dir + "angle.sh", "nalpha=.*", "nalpha=%d" % nalpha)
    # replace (out_dir + "angle.sh", "nangle=.*", "nangle=%d" % nangle)

def gen_rid (out_dir,
              gen_file,
              cv_file,
              mol_dir,
              res_dir,
              rid_dir) :
    mol_dir += "/"
    res_dir += "/"
    rid_dir += "/"
    out_dir += "/"
    conf_file = mol_dir + "conf.gro"
    cv_dim_list = cal_cv_dim (conf_file, cv_file)
    cv_dim = sum(cv_dim_list)
    print ("cv dim:             %d" % cv_dim)

    create_path (out_dir)
    
    # rid root
    for ii in rid_file_copy :        
        if os.path.isdir (rid_dir + ii) :
            shutil.copytree (rid_dir + ii, out_dir + ii)
        elif os.path.isfile (rid_dir + ii) :
            shutil.copy (rid_dir + ii, out_dir)
    assert rid_run in rid_file_copy, "rid file should have run.py"
    assert rid_param in rid_file_copy, "rid file should have param.json"
    assert "template" in rid_file_copy, "rid file should have template"
    # replace (out_dir + rid_run, "cv_dim = .*", "cv_dim = %d" % cv_dim)    
    replace (out_dir + "lib/modeling.py", "cv_dim = .*", "cv_dim = %d" % cv_dim)    
    replace (out_dir + rid_param, 
             "\"template_dir\":.*", "\"template_dir\":\t\"%s\","  % ("./template"))

    # template/00.enhcMD
    replace (out_dir + "/template/00.enhcMD/test.std.py", "cv_dim = .*", "cv_dim = %d" % cv_dim)    
    ret = general_plumed ("dpbias", conf_file, cv_file)
    with open(out_dir + "/template/00.enhcMD/plumed.dat", "w") as fp:
        fp.write (ret)
    ret = general_plumed ("bf", conf_file, cv_file)
    with open(out_dir + "/template/00.enhcMD/plumed.bf.dat", "w") as fp:
        fp.write (ret)

    # template/01.resMD
    gen_res (out_dir + "/template/01.resMD", gen_file, cv_file, mol_dir, res_dir)

    # template/02.train
    replace (out_dir + "/template/02.train/model.py", "cv_dim = .*", "cv_dim = %d" % cv_dim)
    replace (out_dir + "/template/02.train/model.py", "cv_dih_dim = .*", "cv_dih_dim = %d" % cv_dim_list[0])
    replace (out_dir + "/template/02.train/model.py", "cv_dist_dim = .*", "cv_dist_dim = %d" % cv_dim_list[1])

    # template/mol
    gen_bf (out_dir + "/template/mol", gen_file, cv_file, mol_dir)

    # template/tools
    replace (out_dir + "/template/tools/cluster_cv.py", "cv_dih_dim = .*", "cv_dih_dim = %d" % cv_dim_list[0])

def gen_rit (out_dir,
              gen_file,
              cv_file,
              mol_dir,
              res_dir,
              rit_dir) :
    mol_dir += "/"
    res_dir += "/"
    rit_dir += "/"
    out_dir += "/"
    conf_file = mol_dir + "conf.gro"
    cv_dim_list = cal_cv_dim (conf_file, cv_file)
    cv_dim = sum(cv_dim_list)
    print ("cv dim:             %d" % cv_dim)

    create_path (out_dir)

    # rit root
    for ii in rit_file_copy :        
        if os.path.isdir (rit_dir + ii) :
            shutil.copytree (rit_dir + ii, out_dir + ii)
        elif os.path.isfile (rit_dir + ii) :
            shutil.copy (rit_dir + ii, out_dir)
    assert rit_run in rit_file_copy, "rit file should have run.py"
    assert rit_param in rit_file_copy, "rit file should have param.json"
    assert "template" in rit_file_copy, "rit file should have template"
    # replace (out_dir + rit_run, "cv_dim = .*", "cv_dim = %d" % cv_dim)    
    replace (out_dir + "lib/modeling.py", "cv_dim = .*", "cv_dim = %d" % cv_dim)    
    replace (out_dir + rit_param, 
             "\"template_dir\":.*", "\"template_dir\":\t\"%s\","  % ("./template"))

    # template/00.enhcMD
    replace (out_dir + "/template/00.enhcMD/test.std.py", "cv_dim = .*", "cv_dim = %d" % cv_dim)    
    ret = general_plumed ("afed", conf_file, cv_file)
    with open(out_dir + "/template/00.enhcMD/plumed.dat", "w") as fp:
        fp.write (ret)

    # template/01.resMD
    gen_res (out_dir + "/template/01.resMD", gen_file, cv_file, mol_dir, res_dir)

    # template/02.train
    replace (out_dir + "/template/02.train/model.py", "cv_dim = .*", "cv_dim = %d" % cv_dim)
    replace (out_dir + "/template/02.train/model.py", "cv_dih_dim = .*", "cv_dih_dim = %d" % cv_dim_list[0])
    replace (out_dir + "/template/02.train/model.py", "cv_dist_dim = .*", "cv_dist_dim = %d" % cv_dim_list[1])

    # template/mol
    gen_bf (out_dir + "/template/mol", gen_file, cv_file, mol_dir)

    # template/tools
    replace (out_dir + "/template/tools/cluster_cv.py", "cv_dih_dim = .*", "cv_dih_dim = %d" % cv_dim_list[0])

def gen_strm (out_dir,
              json_file,
              mol_dir,
              res_dir,
              strm_dir,
              tools_dir) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)    
    nalpha = jdata['nalpha']
    alpha_idx_fmt = jdata['alpha_idx_fmt']    
    angle_idx_fmt = jdata['angle_idx_fmt']    
    angle_idxes = jdata['angle_idxes']    

    nangle = len(angle_idxes)
    mol_dir += "/"
    res_dir += "/"
    strm_dir += "/"
    tools_dir += "/"
    out_dir += "/"
    cv_dim = nalpha * nangle    

    create_path (out_dir)
    
    angle_idxes_arg=""
    angle_idxes = np.reshape(angle_idxes, [-1])
    for ii in angle_idxes :
        angle_idxes_arg += " " + str(ii)

    # strm root
    for ii in strm_file_copy :        
        if os.path.isdir (strm_dir + ii) :
            shutil.copytree (strm_dir + ii, out_dir + ii)
        elif os.path.isfile (strm_dir + ii) :
            shutil.copy (strm_dir + ii, out_dir)
    replace (out_dir + "param.json", 
             "\"template_dir\":.*", "\"template_dir\":\t\"%s\","  % ("./template"))

    # template/00.resMD
    gen_res (out_dir + "/template/00.resMD", json_file, mol_dir, res_dir, tools_dir)

    # template/mol
    gen_bf (out_dir + "/template/mol", json_file, mol_dir, tools_dir)
    shutil.copy (tools_dir + "conf.ang.sh", out_dir + "/template/mol")
    replace (out_dir + "/template/mol/conf.ang.sh", "nalpha=.*", "nalpha=%d" % nalpha)
    replace (out_dir + "/template/mol/conf.ang.sh", "nangle=.*", "nangle=%d" % nangle)


def _main () :
    parser = argparse.ArgumentParser()
    parser.add_argument("TASK", type=str, 
                        help="the task")
    parser.add_argument("GEN_DEF", type=str, 
                        help="the json database defining the generation")
    parser.add_argument("CV_DEF", type=str, 
                        help="the json database defining CVs")
    parser.add_argument("MOL", type=str, 
                        help="the mol dir")
    parser.add_argument("-o", "--output", type=str, default = "out",
                        help="the output dir")
    args = parser.parse_args()

    # fp = open (args.JSON, 'r')
    # jdata = json.load (fp)    
    # nalpha = jdata["nalpha"]

    base_path = os.path.dirname(os.path.realpath(__file__)) + "/"
    mol_dir = base_path +  args.MOL
    res_dir = base_path + "res"
    afed_dir = base_path + "afed"
    rid_dir = base_path + "abnn"
    rit_dir = base_path + "abnn"
    strm_dir = base_path + "strm"
    tools_dir = base_path + "tools"
    
    print ("task:               %s" % args.TASK)
    print ("using gen def:      %s" % args.GEN_DEF)
    print ("using cv  def:      %s" % args.CV_DEF)
    print ("using mol dir:      %s" % mol_dir)
    print ("output to:          %s" % base_path + args.output)

    if args.TASK == "bf" : 
        gen_bf (args.output, args.GEN_DEF, args.CV_DEF, mol_dir)
    elif args.TASK == "res" :
        gen_res (args.output, args.GEN_DEF, args.CV_DEF, mol_dir, res_dir)
    elif args.TASK == "afed" :
        gen_afed (args.output, args.GEN_DEF, args.CV_DEF, mol_dir, afed_dir)
    elif args.TASK == "abnn" or args.TASK == "rid" :
        gen_rid (args.output, args.GEN_DEF, args.CV_DEF, mol_dir, res_dir, rid_dir)
    elif args.TASK == "rit" :
        gen_rit (args.output, args.GEN_DEF, args.CV_DEF, mol_dir, res_dir, rit_dir)
    # elif args.TASK == "strm" :
    #     gen_strm (args.output, args.JSON, mol_dir, res_dir, strm_dir, tools_dir)

if __name__ == '__main__':
    _main()

