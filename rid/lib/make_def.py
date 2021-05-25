import os, glob
from rid.lib.utils import replace, print_list, print_repeat_list

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

                if not any([(kk in shifted_resid) for kk in atom_name_list]) :
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
                
    return angle_names, angle_atom_idxes

def make_general_dist_def (residues, 
                           residue_atoms,
                           sel_residue_names,
                           sel_atom_names,
                           fmt_residue = "%02d",
                           exclude = 7) :
   
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

    return dist_names, dist_atom_idxes

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
