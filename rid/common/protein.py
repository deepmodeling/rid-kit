import MDAnalysis as mda
from typing import List


def get_all_dihedral_index(file_path: str):
    u = mda.Universe(file_path)
    all_res_list = []
    for seg in u.segments:
        chain_res_list = seg.residues.resindices
        if len(chain_res_list) <= 2:
            continue
        else:
            all_res_list += chain_res_list[1:-1].tolist()
    print("The dihedral angle indexes selected are:", all_res_list)
    return all_res_list


def get_dihedral_info(file_path: str):
    u = mda.Universe(file_path)
    dihedral_angle = {}
    for residue in u.residues:
        dihedral_angle[residue.resid] = {}
        if residue.phi_selection() is not None:
            dihedral_angle[residue.resid]["phi"] = list(residue.phi_selection().ids)
        if residue.psi_selection() is not None:
            dihedral_angle[residue.resid]["psi"] = list(residue.psi_selection().ids)
    return dihedral_angle


def get_dihedral_from_resid(file_path: str, selected_id: List[int]):
    if len(selected_id) == 0:
        return []
    u = mda.Universe(file_path)
    selected_dihedral_angle = {}
    for sid in selected_id:
        residue = u.residues[sid - 1]
        selected_dihedral_angle[sid] = {}
        if residue.phi_selection() is not None:
            selected_dihedral_angle[residue.resid]["phi"] = list(residue.phi_selection().ids)
        if residue.psi_selection() is not None:
            selected_dihedral_angle[residue.resid]["psi"] = list(residue.psi_selection().ids)
          
    return selected_dihedral_angle



