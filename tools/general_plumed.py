#!/usr/bin/env python3

import argparse
import json
import numpy as np
import os, sys
sys.path.append (os.path.dirname(os.path.realpath(__file__)))
from make_ndx import make_ndx
from make_ndx import make_protein_atom_index

def make_general_angle_def (residue_atoms, 
                            dih_angles, 
                            fmt_alpha = "%2d",
                            fmt_angle = "%2d") :
    """
    Inputs:
    residue_atoms:      the atoms in each residule, returned by make_ndx
    dih_angles:         the definition of dihedral angles
    fmt_alpha:          the format of printing residue index
    fmt_angle:          the format of printing angle index
    
    Returns:
    angle_names:        the dihedral angle names in format "resid_idx-angle_idx"
    angle_atom_idxes:   the atom indexs of each dihedral angle
    """
    angle_names = []
    angle_atom_idxes = []
    for ii in range(len(residue_atoms)) :
        resid = residue_atoms[ii]
        for jj in range(len(dih_angles)) :
            angle = dih_angles[jj]            
            angle_print = "dih-" + (fmt_alpha % ii) + "-" + (fmt_angle % jj)
            find_angle = True
            atom_idxes = []
            for atom in angle :
                shifted_resid_idx = ii + atom['resid_shift']
                if shifted_resid_idx < 0 or shifted_resid_idx == len(residue_atoms) :
                    find_angle = False 
                    break
                shifted_resid = residue_atoms[shifted_resid_idx]
                atom_name_list = atom['name']
                if not any([(ii in shifted_resid) for ii in atom_name_list]) :
                    find_angle = False
                    break
                for atom_name in atom_name_list :
                    if atom_name in shifted_resid :
                        atom_idxes.append(shifted_resid[atom_name])
                        break
            if find_angle :
                assert(len(atom_idxes) == 4)
                angle_names.append(angle_print)
                angle_atom_idxes.append(atom_idxes)
                # mylist = ""
                # for kk in atom_idxes:
                #     if len(mylist) == 0 :
                #         mylist = str(kk)
                #     else :
                #         mylist += "," + str(kk)
                # ret += ( (angle_print + ":") + " " + 
                #          "TORSION " + 
                #          "ATOMS=" + 
                #          mylist + " " +
                #          "\n")
    return angle_names, angle_atom_idxes

def make_general_dist_def (residues, 
                           residue_atoms,
                           sel_residue_names,
                           sel_atom_names,
                           fmt_residue = "%02d",
                           exclude = 7) :
    # print (residues)
    # print (sel_residue_names)    
    sel_residue_idx = []
    sel_atom_idx = []
    for ii in range(len(residues)) :
        if residues[ii][0] in sel_residue_names : 
            find_atom = False
            for jj in sel_atom_names : 
                if jj in residue_atoms[ii] :
                    find_atom = True
                    break
            if not find_atom : 
                continue
            sel_atom_name = jj
            # print (ii, residues[ii], residue_atoms[ii], residue_atoms[ii][sel_atom_name])
            sel_residue_idx.append(ii)
            sel_atom_idx.append(residue_atoms[ii][sel_atom_name])

    dist_names = []
    dist_atom_idxes = []
    for ii in range(len(sel_residue_idx)) :
        for jj in range(ii+1, len(sel_residue_idx)) :
            ri = sel_residue_idx[ii]
            rj = sel_residue_idx[jj]
            ai = sel_atom_idx[ii]
            aj = sel_atom_idx[jj]
            if rj - ri < exclude : 
                continue
            dist_names.append ("dist-" + (fmt_residue % ri) + "-" + (fmt_residue % rj))
            dist_atom_idxes.append ([ai, aj])
    # print (dist_names)
    # print (dist_atom_idxes)
    return dist_names, dist_atom_idxes

def print_list (tmp, 
                suffix = "") :
    mylist = ""
    for kk in tmp:
        if len(mylist) == 0 :
            mylist = str(kk) + suffix
        else :
            mylist += "," + str(kk) + suffix
    return mylist

def print_repeat_list (numb, item) :
    mylist=""
    for ii in range(numb) :
        if ii == 0: 
            mylist = str(item)
        else :
            mylist += "," + str(item)
    return mylist

def make_angle_def (angle_names,
                    angle_atom_idxes) :
    ret = ""
    for angle_print, atom_idxes in zip(angle_names, angle_atom_idxes) :
        mylist = print_list(atom_idxes)
        ret += ( (angle_print + ":") + " " + 
                "TORSION " + 
                "ATOMS=" + 
                mylist + " " +
                "\n")    
    return ret

def make_dist_def (dist_names, 
                   dist_atom_idxes) :
    ret = ""
    for dist_print, atom_idxes in zip(dist_names, dist_atom_idxes) :
        mylist = print_list(atom_idxes)
        ret += (dist_print + ":" + " " + 
                "DISTANCE" + " "
                "ATOMS=" + 
                mylist + " " +
                "\n")    
    return ret

def make_restraint (angle_names, 
                    kappa,
                    at) :
    ret = ""
    res_names = []
    for angle_print in angle_names :
        res_name = "res-" + angle_print
        res_names.append(res_name)
        ret += ( (res_name + ":") + " " + 
                 "RESTRAINT " + 
                 ("ARG=" + angle_print).ljust(16) + " " + 
                 "KAPPA=" + str(kappa) + " " +
                 "AT=" + str(at) + " " + 
                 "\n" )
    return ret, res_names

def make_afed (angle_names,
               temp,
               kappa,
               tau,
               gamma) :
    arg_list = print_list (angle_names)
    arg_list_fict = print_list (angle_names, "_fict")
    nargs = len(angle_names)
    kappa_list = print_repeat_list (nargs, kappa)
    gamma_list = print_repeat_list (nargs, gamma)
    tau_list = print_repeat_list (nargs, tau)
    return ("ex: " + 
            "EXTENDED_LAGRANGIAN " + 
            "TEMP=" + str(temp) + " " +
            "ARG=" + arg_list + " " +
            "KAPPA=" + kappa_list + " " + 
            "TAU=" + tau_list + " " + 
            "FRICTION=" + gamma_list + " " + 
            "\n")

def make_deep_bias (angle_names,
                    trust_lvl_1 = 1.0,
                    trust_lvl_2 = 2.0,
                    model = "graph.pb") :
    arg_list = print_list (angle_names)
    return ("dpfe: " + 
            "DEEPFE " + 
            "TRUST_LVL_1=" + str(trust_lvl_1) + " " + 
            "TRUST_LVL_2=" + str(trust_lvl_2) + " " +
            "MODEL=" + model + " " +
            "ARG=" + str(arg_list) + " " + 
            "\n" )

def make_print (names,
                stride,
                file_name) :
    arg_list = print_list(names)
    return ("PRINT " +
            "STRIDE=" + str(stride) + " " +
            "ARG=" + str(arg_list) + " " + 
            "FILE=" + str(file_name) + " " + 
            "\n" )

def make_wholemolecules (atom_index) :
    arg_list = print_list(atom_index)
    return ("WHOLEMOLECULES" + " " +
            "ENTITY0=" + arg_list +
            "\n")            

def cal_cv_dim (conf_file, cv_file) :    
    cfile = conf_file
    jfile = cv_file
    residues, residue_atoms = make_ndx (cfile)
    fp = open (jfile, 'r')
    jdata = json.load (fp)
    dih_angles = jdata["dih_angles"]
    fmt_alpha = jdata["alpha_idx_fmt"]
    fmt_angle = jdata["angle_idx_fmt"]
    hp_residues = []
    dist_atom = []
    dist_excl = 10000
    if "hp_residues" in jdata :
        hp_residues = jdata["hp_residues"]
    if "dist_atom" in jdata:
        dist_atom = jdata["dist_atom"]
    if "dist_excl" in jdata:
        dist_excl = jdata["dist_excl"]

    angle_names, angle_atom_idxes = make_general_angle_def(residue_atoms, dih_angles, fmt_alpha, fmt_angle)
    dist_names, dist_atom_idxes = make_general_dist_def(residues, residue_atoms, hp_residues, dist_atom, fmt_alpha, dist_excl)
    return [len(angle_names), len(dist_names)]

def general_plumed (TASK,
                    CONF,
                    JSON,
                    kappa = 500.0,
                    temp = 3000.0,
                    tau = 10.0,
                    gamma = 0.1,
                    pstride = 5,
                    pfile = "plm.out") :
    residues, residue_atoms = make_ndx (CONF)
    protein_atom_idxes = make_protein_atom_index(CONF)
    fp = open (JSON, 'r')
    jdata = json.load (fp)
    dih_angles = jdata["dih_angles"]
    fmt_alpha = jdata["alpha_idx_fmt"]
    fmt_angle = jdata["angle_idx_fmt"]
    hp_residues = []
    dist_atom = []
    dist_excl = 10000
    if "hp_residues" in jdata :
        hp_residues = jdata["hp_residues"]
    if "dist_atom" in jdata:
        dist_atom = jdata["dist_atom"]
    if "dist_excl" in jdata:
        dist_excl = jdata["dist_excl"]

    angle_names, angle_atom_idxes = make_general_angle_def(residue_atoms, dih_angles, fmt_alpha, fmt_angle)

    dist_names, dist_atom_idxes = make_general_dist_def(residues, residue_atoms, hp_residues, dist_atom, fmt_alpha, dist_excl)

    ret = ""
    if len(dist_names) > 0 :
        ret += make_wholemolecules(protein_atom_idxes)
        ret += "\n"
    ret += (make_angle_def(angle_names, angle_atom_idxes))
    ret += (make_dist_def(dist_names, dist_atom_idxes))
    ret += "\n"

    cv_names = angle_names + dist_names
    if TASK == "res" :
        ptr, ptr_names = make_restraint(cv_names, kappa, 0.0)
        ret += (ptr)
        ret += "\n"
    elif TASK == "afed" :
        ret += (make_afed(cv_names, temp, kappa, tau, gamma))
        ret += "\n"
    elif TASK == "dpbias" :
        ret += (make_deep_bias(cv_names))
        ret += "\n"
    elif TASK == "bf" :
        None
    else :
        raise RuntimeError ("unknow task: " + TASK)
    ret += (make_print(cv_names, pstride, pfile))    
    ret += "\n"
    return ret

def _main () :
    parser = argparse.ArgumentParser()
    parser.add_argument("TASK", type = str,
                        help="the type of task, either res, afed or dpbias")
    parser.add_argument("CONF", type=str, 
                        help="the conf file")
    parser.add_argument("JSON", type=str, 
                        help="the json file defining the dih angles")
    parser.add_argument("-k", "--kappa", default=500, type=float, 
                        help="the spring constant")
    parser.add_argument("-T", "--temp", default=3000.0, type=float, 
                        help="the temperature of afed")
    parser.add_argument("-t", "--tau", default=10.0, type=float, 
                        help="the relaxation timescale of afed")
    parser.add_argument("-g", "--gamma", default=0.1, type=float, 
                        help="the frection const of afed")
    parser.add_argument("-s", "--stride", default=5, type=int, 
                        help="the printing stride")
    parser.add_argument("-f", "--print-file", default="plm.out", type=str, 
                        help="the printing file")
    args = parser.parse_args()

    ret = general_plumed (args.TASK, args.CONF, args.JSON, args.kappa, args.temp, args.tau, args.gamma, args.stride, args.print_file)

    print (ret)

if __name__ == '__main__':
    _main()
    
