import os
import sys
from typing import List, Dict
import logging
import MDAnalysis as mda


logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

def get_all_dihedral_index(file_path: str):
    u = mda.Universe(file_path)
    all_res_list = []
    for seg in u.segments:
        chain_res_list = seg.residues.resindices
        if len(chain_res_list) <= 2:
            continue
        else:
            all_res_list += chain_res_list[1:-1].tolist()
    logger.debug("The dihedral angle indexes selected are:", all_res_list)
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


def get_dihedral_from_resid(file_path: str, selected_resid: List[int]) -> Dict:
    if len(selected_resid) == 0:
        return {}
    u = mda.Universe(file_path)
    selected_dihedral_angle = {}
    for sid in selected_resid:
        residue = u.residues[sid - 1]
        selected_dihedral_angle[sid] = {}
        if residue.phi_selection() is not None:
            selected_dihedral_angle[residue.resid]["phi"] = list(residue.phi_selection().ids)
        if residue.psi_selection() is not None:
            selected_dihedral_angle[residue.resid]["psi"] = list(residue.psi_selection().ids)
    num_cv = len(selected_dihedral_angle.keys())
    logger.info(f"{num_cv} CVs have been created.")
    return selected_dihedral_angle



