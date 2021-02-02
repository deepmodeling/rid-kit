#!/usr/bin/env python3

import argparse
import numpy as np

def make_angle_def (idx_cum, 
                    angle_idxes,
                    fmt_alpha,
                    fmt_angle) :
    ret = ""
    for ii in range(len(idx_cum)) :
        for jj in range(len(angle_idxes)) :
            angle_print = (fmt_alpha % ii) + "-" + (fmt_angle % jj)
            mylist = ""
            for kk in angle_idxes[jj] :
                if len(mylist) == 0 :
                    mylist = str(idx_cum[ii] + kk)
                else :
                    mylist += "," + str(idx_cum[ii] + kk)
            ret += ( (angle_print + ":") + " " + 
                     "TORSION " + 
                     "ATOMS=" + 
                     mylist + " " +
                     "\n")
    return ret

def make_arg_list (idx_cum, 
                   angle_idxes,
                   fmt_alpha,
                   fmt_angle, 
                   suffix = "") :
    arg_list=""
    for ii in range(len(idx_cum)) :
        for jj in range(len(angle_idxes)) :
            angle_print = (fmt_alpha % ii) + "-" + (fmt_angle % jj)
            if len(arg_list) == 0: 
                arg_list = angle_print + suffix
            else :
                arg_list += "," + angle_print + suffix
    return arg_list

def make_restraint (idx_cum, 
                    angle_idxes,
                    fmt_alpha,
                    fmt_angle,
                    kappa,
                    at) :
    ret = ""
    for ii in range(len(idx_cum)) :
        for jj in range(len(angle_idxes)) :
            angle_print = (fmt_alpha % ii) + "-" + (fmt_angle % jj)
            ret += ( ("res-" + angle_print + ":") + " " + 
                     "RESTRAINT " + 
                     ("ARG=" + angle_print).ljust(16) + " " + 
                     "KAPPA=" + str(kappa) + " " +
                     "AT=" + str(at) + " " + 
                     "\n" )
    return ret

def make_repeat_list (numb, item) :
    mylist=""
    for ii in range(numb) :
        if ii == 0: 
            mylist = str(item)
        else :
            mylist += "," + str(item)
    return mylist

def make_afed (idx_cum,
               angle_idxes,
               fmt_alpha,
               fmt_angle,
               temp,
               kappa,
               tau,
               gamma) :
    arg_list = make_arg_list (idx_cum, angle_idxes, fmt_alpha, fmt_angle)
    arg_list_fict = make_arg_list (idx_cum, angle_idxes, fmt_alpha, fmt_angle, "_fict")
    nargs = len(idx_cum) * len(angle_idxes)
    kappa_list = make_repeat_list (nargs, kappa)
    gamma_list = make_repeat_list (nargs, gamma)
    tau_list = make_repeat_list (nargs, tau)
    return ("ex: " + 
            "EXTENDED_LAGRANGIAN " + 
            "TEMP=" + str(temp) + " " +
            "ARG=" + arg_list + " " +
            "KAPPA=" + kappa_list + " " + 
            "TAU=" + tau_list + " " + 
            "FRICTION=" + gamma_list + " ")

def make_deep_bias (idx_cum,
                    angle_idxes,
                    fmt_alpha,
                    fmt_angle,
                    trust_lvl_1 = 1.0,
                    trust_lvl_2 = 2.0,
                    model = "graph.pb") :
    arg_list = make_arg_list (idx_cum, angle_idxes, fmt_alpha, fmt_angle)
    return ("dpfe: " + 
            "DEEPFE " + 
            "TRUST_LVL_1=" + str(trust_lvl_1) + " " + 
            "TRUST_LVL_2=" + str(trust_lvl_2) + " " +
            "MODEL=" + model + " " +
            "ARG=" + str(arg_list) + " " )

def make_print (idx_cum, 
                angle_idxes,
                fmt_alpha,
                fmt_angle,
                stride,
                file_name) :
    arg_list = make_arg_list (idx_cum, angle_idxes, fmt_alpha, fmt_angle)
    return ("PRINT " +
            "STRIDE=" + str(stride) + " " +
            "ARG=" + str(arg_list) + " " + 
            "FILE=" + str(file_name) + " " )

def make_shift_cum (numb_calpha) :
    idx_shift = np.ones(numb_calpha, dtype=np.int) * 10
    idx_shift[0] = 0
    idx_cum = np.cumsum (idx_shift)
    return idx_cum

def _main () :
    parser = argparse.ArgumentParser()
    parser.add_argument("TASK", type = str,
                        help="the type of task, either res, afed or dpbias")
    parser.add_argument("-a", "--numb-alpha", default=0, type=int, 
                        help="the number of alpha carbons")
    parser.add_argument("-i", "--angle-indexes", default=[], nargs= '*', type=int, 
                        help="the indexes of dihedral angle atoms")
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
    parser.add_argument("--alpha-fmt", default="%02d", type=str, 
                        help="format of alpha")
    parser.add_argument("--angle-fmt", default="%02d", type=str, 
                        help="format of angle")

    args = parser.parse_args()
    
    numb_calpha = args.numb_alpha
    # shape should be numb_angles x 4
    angle_idxes = np.reshape (np.array(args.angle_indexes, dtype=np.int), [-1,4])    

    print ("# number of alpha carbon(s): " + str(numb_calpha))
    if numb_calpha == 0: 
        return

    idx_cum = make_shift_cum (numb_calpha)
    
    print (make_angle_def      (idx_cum, angle_idxes, args.alpha_fmt, args.angle_fmt))
    if args.TASK == "res" :
        print (make_restraint  (idx_cum, angle_idxes, args.alpha_fmt, args.angle_fmt, args.kappa, 0.0))
        print ("")
    elif args.TASK == "afed" :
        print (make_afed       (idx_cum, angle_idxes, args.alpha_fmt, args.angle_fmt, args.temp, args.kappa, args.tau, args.gamma))
        print ("")
    elif args.TASK == "dpbias" :
        print (make_deep_bias  (idx_cum, angle_idxes, args.alpha_fmt, args.angle_fmt))
        print ("")
    else :
        raise RuntimeError ("unknow task: " + args.TASK)
    print (make_print          (idx_cum, angle_idxes, args.alpha_fmt, args.angle_fmt, args.stride, args.print_file))

if __name__ == '__main__':
    _main()
    
