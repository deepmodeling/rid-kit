from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Parameter
)

from typing import List, Optional, Union, Dict
from pathlib import Path
import numpy as np
from rid.constants import (
        cv_force_out
    )
from rid.common.gromacs.trjconv import generate_coords, generate_forces
from rid.utils import load_txt, set_directory
from rid.constants import gmx_coord_name, gmx_force_name,f_cvt,kb,mf_fig
from matplotlib import pyplot as plt
import os
from parmed import gromacs
import parmed as pmd

def plot_mf_average_all(mf_average, task_path, dt):
    plt.figure(figsize=(10, 8), dpi=100)
    xlist = [i*dt for i in range(len(mf_average))]
    for index in range(mf_average.shape[1]):
        plt.scatter(xlist, mf_average[:,index],label="d%s"%(index+1))
    plt.title("running average of mean force")
    plt.xlabel("constrained simulation time (ps)")
    plt.ylabel("mean force (KJ/(mol*nm))")
    plt.legend()
    plt.savefig(task_path.joinpath(mf_fig))

def phase_factor(r_cv, cv_init, selected_atomid_simple, mass_list):
    C_m = np.zeros(shape=(len(cv_init), len(cv_init)))
    C_list = np.zeros(shape=(len(cv_init), len(r_cv)))
    for index in range(len(cv_init)):
        atom_id1 = selected_atomid_simple[index][0]
        atom_id2 = selected_atomid_simple[index][1]
        C_list[index][atom_id1*3] = (r_cv[atom_id1*3] - r_cv[atom_id2*3])/cv_init[index]
        C_list[index][atom_id1*3+1] = (r_cv[atom_id1*3+1] - r_cv[atom_id2*3+1])/cv_init[index]
        C_list[index][atom_id1*3+2] = (r_cv[atom_id1*3+2] - r_cv[atom_id2*3+2])/cv_init[index]
        C_list[index][atom_id2*3] = -(r_cv[atom_id1*3] - r_cv[atom_id2*3])/cv_init[index]
        C_list[index][atom_id2*3+1] = -(r_cv[atom_id1*3+1] - r_cv[atom_id2*3+1])/cv_init[index]
        C_list[index][atom_id2*3+2] = -(r_cv[atom_id1*3+2] - r_cv[atom_id2*3+2])/cv_init[index]
    
    Mass = np.diag(mass_list)
    Mass = np.linalg.inv(Mass)
    for index_i in range(len(cv_init)):
        for index_j in range(len(cv_init)):
            C_m[index_i,index_j] = np.dot(np.dot(C_list[index_i],Mass),C_list[index_j])
    
    A = np.linalg.det(C_m)
    return A

def pseudo_inv(r_cv, cv_init, selected_atomid_simple):
    A = np.zeros((len(r_cv),len(cv_init)))
    for index in range(len(cv_init)):
        atom_id1 = selected_atomid_simple[index][0]
        atom_id2 = selected_atomid_simple[index][1]
        A[atom_id1*3,index] = (r_cv[atom_id1*3] - r_cv[atom_id2*3])/cv_init[index]
        A[atom_id1*3+1,index] = (r_cv[atom_id1*3+1] - r_cv[atom_id2*3+1])/cv_init[index]
        A[atom_id1*3+2,index] = (r_cv[atom_id1*3+2] - r_cv[atom_id2*3+2])/cv_init[index]
        A[atom_id2*3,index] = -(r_cv[atom_id1*3] - r_cv[atom_id2*3])/cv_init[index]
        A[atom_id2*3+1,index] = -(r_cv[atom_id1*3+1] - r_cv[atom_id2*3+1])/cv_init[index]
        A[atom_id2*3+2,index] = -(r_cv[atom_id1*3+2] - r_cv[atom_id2*3+2])/cv_init[index]
    
    U, S, Vh = np.linalg.svd(A, full_matrices = False)
    B = np.matmul(np.matmul(np.transpose(Vh),np.linalg.inv(np.diag(S))),np.transpose(U))
    
    return B

if os.path.exists("/gromacs/share/gromacs/top"):
    gromacs.GROMACS_TOPDIR = "/gromacs/share/gromacs/top"
elif os.path.exists("/opt/conda/share/gromacs/top"):
    gromacs.GROMACS_TOPDIR = "/opt/conda/share/gromacs/top"

def CalcMF(
    conf: str,
    task_name: str,
    plm_out: str,
    cv_config: dict,
    label_config: dict,
    tail: float = 0.9,
    topology: Optional[str] = None,
    frame_coords: Optional[str] = None,
    frame_forces: Optional[str] = None,
    at: Optional[str] = None
):
    """
    Calculate mean force with the results of restrained or constrained MD. 

    CalcMF will handle periodic CVs by `angular_mask`. 
    To get the mean value of CVs near equilibrium state under restrained or constrained MD, only part of outputs of CV values 
    (the last `tail` part, for example, the last 90% CV values) are used. 
    """
    
    if label_config["method"] == "restrained":
        data = load_txt(plm_out)
        data = data[:, 1:]  # removr the first column(time index).
        centers = data[0,:]

        nframes = data.shape[0]
        
        angular_boolean = (np.array(cv_config["angular_mask"], dtype=int) == 1)
        init_angle = data[0][angular_boolean]
        for ii in range(1, nframes):
            current_angle = data[ii][angular_boolean]
            angular_diff = current_angle - init_angle
            current_angle[angular_diff < -np.pi] += 2 * np.pi
            current_angle[angular_diff >= np.pi] -= 2 * np.pi
            data[ii][angular_boolean] = current_angle

        start_f = int(nframes * (1-tail))
        avgins = np.average(data[start_f:, :], axis=0)
        
        mf_avg_list = []
        for simu_frames in range(int(nframes*0.1), nframes,1):
            start_f = int(simu_frames*(1-0.9))
            avgins_ = np.average(data[start_f:simu_frames, :], axis=0)
            diff_ = avgins_ - centers
            angular_diff_ = diff_[angular_boolean]
            angular_diff_[angular_diff_ < -np.pi] += 2 * np.pi
            angular_diff_[angular_diff_ >  np.pi] -= 2 * np.pi
            diff_[angular_boolean] = angular_diff_
            ff_ = np.multiply(label_config["kappas"], diff_)
            mf_avg_list.append(ff_)

        mf_avg_list = np.array(mf_avg_list)
        mid_f = int(nframes * 0.5)
        mf_std = np.std(mf_avg_list[mid_f:,:], axis=0)
        
        diff = avgins - centers
        angular_diff = diff[angular_boolean]
        angular_diff[angular_diff < -np.pi] += 2 * np.pi
        angular_diff[angular_diff >  np.pi] -= 2 * np.pi
        diff[angular_boolean] = angular_diff
        ff = np.multiply(label_config["kappas"], diff)
        
        task_path = Path(task_name)
        task_path.mkdir(exist_ok=True, parents=True)
        cv_forces = np.concatenate((centers, ff))
        np.savetxt(task_path.joinpath(cv_force_out),  np.reshape(cv_forces, [1, -1]), fmt='%.10f')
        plot_mf_average_all(mf_avg_list, task_path, dt = label_config["dt"]*label_config["output_freq"])
        
        with open(task_path.joinpath("mf_info.out"),"w") as f:
            f.write("cv list value      ")
            for cv_ in centers:
                f.write("%.4f "%cv_)
            f.write("\n")
            f.write("mean force value   ")
            for mf_ in ff:
                f.write("%.4f "%mf_)
            f.write("\n")
            f.write("mean force std     ")
            for mf_std_ in mf_std:
                f.write("%.4f "%mf_std_)
                
    elif label_config["method"] == "constrained":
        data = load_txt(plm_out)
        data = data[:, 1:]  # removr the first column(time index).
        centers = data[0,:]
        
        # Kb to KJ/mol
        KB = kb*f_cvt
        T = float(label_config["temperature"])
        coords = np.loadtxt(frame_coords,comments=["#","@"])
        
        if "units" in cv_config:
            length_units = cv_config["units"]
        else:
            length_units = None
        if length_units == None or length_units == "nm":
            # coords units nm
            coords = coords[:,1:]
            forces = np.loadtxt(frame_forces,comments=["#","@"])
            # forces units to KJ/(mol*nm)
            forces = forces[:,1:]
        elif length_units == "A" or length_units == "Angstrom":
            # coords units nm
            coords = coords[:,1:]*10
            forces = np.loadtxt(frame_forces,comments=["#","@"])
            # forces units to KJ/(mol*nm)
            forces = forces[:,1:]/10
        else:
            raise ValueError("Valid length units must be nm or A")
        
        mflist = []
        mflist_phase = []
        phase_list = []
        mf_average_phase = []
        selected_atomid = cv_config["selected_atomid"]
        selected_atoms = list(set([item for sublist in selected_atomid for item in sublist]))
        selected_atomid_simple = []
        for atom_pairs in selected_atomid:
            selected_atomid_simple.append([selected_atoms.index(i) for i in atom_pairs])
            
        # calculate mass matrix of the system
        system = pmd.load_file(os.path.abspath(topology))
        mass_list = [system.atoms[i].mass for i in range(len(system.atoms))]
        mass_list_simple = []
        for atom_id in selected_atoms:
            atom_id -= 1
            mass_list_simple.append(mass_list[atom_id])
            mass_list_simple.append(mass_list[atom_id])
            mass_list_simple.append(mass_list[atom_id])
        
        # calculate mean force via singular value decomposition(SVD)
        eps = 0.0001
        for index in range(np.shape(coords)[0]):
            coordlist = coords[index]
            forcelist = forces[index]
            r_cv = []
            f_cv = []
            for atom_id in selected_atoms:
                atom_id -= 1
                r_cv.append(coordlist[atom_id*3])
                r_cv.append(coordlist[atom_id*3+1])
                r_cv.append(coordlist[atom_id*3+2])
                f_cv.append(forcelist[atom_id*3])
                f_cv.append(forcelist[atom_id*3+1])
                f_cv.append(forcelist[atom_id*3+2])
            B = pseudo_inv(r_cv, centers, selected_atomid_simple)
            # print(B.shape)
            dBdx = []
            for index in range(len(r_cv)):
                r1 = r_cv.copy()
                r1[index] += eps
                B1 = pseudo_inv(r1,centers,selected_atomid_simple)
                r2 = r_cv.copy()
                r2[index] -= eps
                B2 = pseudo_inv(r2, centers,selected_atomid_simple)
                dBdx.append((B1-B2)/(2*eps))
            dBdx = np.array(dBdx)
            phase = phase_factor(r_cv, centers, selected_atomid_simple, mass_list_simple)
            mf = np.matmul(B,f_cv) + KB*T*np.trace(dBdx, axis1=0, axis2=2)
            
            phase_list.append(1/(phase**0.5))
            mflist.append(mf)
            mflist_phase.append(mf/(phase**0.5))
        
        mflist = np.array(mflist)
        mflist_phase = np.array(mflist_phase)
        nframes = np.shape(coords)[0]            
        
        start_f = int(nframes * (1-tail))
        mid_f = int(nframes * 0.5)
        mf_avg_without_norm = np.average(mflist_phase[start_f:, :], axis=0)
        phase_avg = np.average(phase_list[start_f:])
        mf_avg_without_phase = np.average(mflist[start_f:, :], axis=0)
        avg_force = mf_avg_without_norm/phase_avg
        
        for simu_frames in range(int(nframes*0.1), nframes,1):
            start_f = int(simu_frames*(1-0.9))
            mf_average_phase.append(np.average(mflist_phase[start_f:simu_frames, :], axis=0)/np.average(phase_list[start_f:simu_frames]))
        
        # calculate mean force statistics
        mf_average_phase = np.array(mf_average_phase)
        mf_difference = avg_force - mf_avg_without_phase
        mf_std = np.std(mf_average_phase[mid_f:,:], axis=0)
        
        task_path = Path(task_name)
        task_path.mkdir(exist_ok=True, parents=True)
        cv_forces = np.concatenate((centers, avg_force))
        np.savetxt(task_path.joinpath(cv_force_out),  np.reshape(cv_forces, [1, -1]), fmt='%.10f')
        plot_mf_average_all(mf_average_phase, task_path, dt = label_config["dt"]*label_config["output_freq"])
        
        with open(task_path.joinpath("mf_info.out"),"w") as f:
            f.write("cv list value      ")
            for cv_ in centers:
                f.write("%.4f "%cv_)
            f.write("\n")
            f.write("mean force value   ")
            for mf_ in avg_force:
                f.write("%.4f "%mf_)
            f.write("\n")
            f.write("diff with phase    ")
            for mf_diff_ in mf_difference:
                f.write("%.4f "%mf_diff_)
            f.write("\n")
            f.write("mean force std     ")
            for mf_std_ in mf_std:
                f.write("%.4f "%mf_std_)
        
    mf_info = None
    if os.path.exists(task_path.joinpath("mf_info.out")):
        mf_info = task_path.joinpath("mf_info.out")
            
    op_out = {
                "cv_forces": task_path.joinpath(cv_force_out),
                "mf_info": mf_info
            }

    return op_out