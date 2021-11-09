import json
from rid.lib.make_ndx import make_ndx
from rid.lib.make_def import make_general_angle_def, make_general_dist_def


def cal_cv_dim_json(conf_file, cv_file):
    cfile = conf_file
    jfile = cv_file
    residues, residue_atoms = make_ndx(cfile)
    fp = open(jfile, 'r')
    jdata = json.load(fp)
    fp.close()
    dih_angles = jdata["dih_angles"]
    fmt_alpha = jdata["alpha_idx_fmt"]
    fmt_angle = jdata["angle_idx_fmt"]
    hp_residues = []
    dist_atom = []
    dist_excl = 10000
    if "hp_residues" in jdata:
        hp_residues = jdata["hp_residues"]
    if "dist_atom" in jdata:
        dist_atom = jdata["dist_atom"]
    if "dist_excl" in jdata:
        dist_excl = jdata["dist_excl"]
    angle_names, angle_atom_idxes = make_general_angle_def(
        residue_atoms, dih_angles, fmt_alpha, fmt_angle)
    dist_names, dist_atom_idxes = make_general_dist_def(
        residues, residue_atoms, hp_residues, dist_atom, fmt_alpha, dist_excl)
    if "selected_index" in jdata:
        selected_index = jdata["selected_index"]
        selected_angle_index = []
        for ssi in selected_index:
            if ssi == 0:
                selected_angle_index.append(0)
            elif ssi != 0:
                selected_angle_index.append(ssi*2-1)
                selected_angle_index.append(ssi*2)

        newselected_angle_index = [
            i for i in selected_angle_index if i < len(angle_names)]
        new_angle_names = [angle_names[i] for i in newselected_angle_index]
        new_angle_atom_idxes = [angle_atom_idxes[i]
                                for i in newselected_angle_index]
    else:
        new_angle_names = angle_names
        new_angle_atom_idxes = angle_atom_idxes

    return [len(new_angle_names), len(dist_names)]

def cal_cv_dim_user(conf_file, cv_file):
    dih_angle_counts = 0
    non_angle_counts = 0
    with open(cv_file, 'r') as fp:
        for line in fp.readlines():
            if ("PRINT" in line) and ("#" not in line):
                print_content = line + "\n"
                break
            if (":" in line) and ("#" not in line):
                if line.split()[1] == "TORSION":
                    dih_angle_counts += 1
                else:
                    non_angle_counts += 1
            elif ("LABEL" in line) and ("#" not in line):
                non_angle_counts += 1
    return [int(dih_angle_counts), int(non_angle_counts)]


def cal_cv_dim(conf_file, cv_file):
    if cv_file.split(".")[-1] == "josn":
        return cal_cv_dim_json(conf_file, cv_file)
    else:
        return cal_cv_dim_user(conf_file, cv_file)