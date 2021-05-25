#!/usr/bin/env python3

import json

amoni_acids = ["ARG", "HIS", "LYS", 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'SEC', 'GLY', 'PRO', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP', 'ACE', 'NME']
main_chain_atom_name = ["N", "CA", "C", "O", "OC1", "CB", "CH3"]

def get_res_idx (line) :
    res = line[0:5]
    res = res.strip()
    return int(res)

def get_res_name (line) :
    res = line[5:10]
    res = res.strip()
    return str(res)

def get_atom_idx (line) :
    atom = line[15:20]
    atom = atom.strip()
    return int(atom)

def get_atom_name (line) :
    atom = line[10:15]
    atom = atom.strip()
    return str(atom)

def make_residue_atoms (res, istart, iend) :
    atom_names = []
    atom_idxes = []
    for ii in range(istart, iend) :
        if get_atom_name(res[ii]) in main_chain_atom_name :
            atom_names.append(get_atom_name(res[ii]))
            atom_idxes.append(get_atom_idx (res[ii]))
    if len(set(atom_names)) != len(atom_names) :
        raise RuntimeError("find duplicated atoms in residue with atoms: %d - %d " % (istart, iend))
    resid = {}
    for nn,ii in zip(atom_names, atom_idxes) :
        resid[nn] = ii
    return resid

def make_ndx (fname) :
    # get lines
    with open(fname) as f:
        content = f.readlines()
    # get atom lines
    res = []
    for ii in range(2,len(content)-1) : 
        res.append(content[ii])
    # check number of residue
    resid_idx = [get_res_idx(ii) for ii in res]
    resid_idx.sort()
    numb_resid = resid_idx[-1]
    # residues atom count
    resid_atom_count = []
    for ii in range(1, numb_resid+1) :
        resid_atom_count.append( sum([idx == ii for idx in resid_idx]) )    
    # construct residues
    residues = []
    residue_atoms = []
    istart = 0
    for numb_atoms in resid_atom_count :
        iend = istart + numb_atoms
        residues.append([get_res_name(res[istart]), istart, iend])
        residue_atoms.append(make_residue_atoms(res, istart, iend))
        istart = iend
    # 
    return residues, residue_atoms

def make_protein_atom_index (fname) :
    ret = []
    with open(fname) as f:
        content = f.readlines()
    res = []
    for ii in range(2,len(content)-1) : 
        res.append(content[ii])
    for ii in res :
        if get_res_name(ii) in amoni_acids :
            ret.append(get_atom_idx(ii))
    return ret

def _main () :
    residues, residue_atoms = make_ndx ("conf.gro")
    fp = open ("dih.json", 'r')
    jdata = json.load (fp)
    dih_angles = jdata["dih_angles"]
    fmt_alpha = jdata["alpha_idx_fmt"]
    fmt_angle = jdata["angle_idx_fmt"]
    
    angle_names, angle_atom_idxes = make_general_angle_def(residue_atoms, dih_angles, fmt_alpha, fmt_angle)
    for angle_print, atom_idxes in zip(angle_names, angle_atom_idxes) :
        mylist = ""
        for kk in atom_idxes:
            if len(mylist) == 0 :
                mylist = str(kk)
            else :
                mylist += "," + str(kk)
        print ( (angle_print + ":") + " " + 
                "TORSION " + 
                "ATOMS=" + 
                mylist + " " +
                "")    

if __name__ == '__main__':
    _main()


