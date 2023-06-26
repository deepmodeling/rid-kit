import os
import sys
from typing import List, Dict, Sequence, Union
import logging
import mdtraj as md
from mdtraj.geometry.dihedral import _atom_sequence, PHI_ATOMS, PSI_ATOMS
import numpy as np
from rid.constants import sel_gro_name_gmx, sel_gro_name
from rid.common.gromacs.trjconv import slice_trjconv

logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def _zip_dict(resi_indices, atom_indices):
    dihedral_dict = {}
    for idx, resi_idx in enumerate(resi_indices):
        dihedral_dict[resi_idx] = atom_indices[idx].tolist()
    return dihedral_dict

def distance(r1,r2):
    d = 0
    for index in range(3):
        d += (r1[index] - r2[index])**2
    d = np.sqrt(d)
    return d

def get_dihedral_info(file_path: str):
    traj = md.load(file_path)
    top = traj.topology
    psi_found_indices, psi_atom_indices = _atom_sequence(top, PSI_ATOMS)
    psi_info = _zip_dict(psi_found_indices + 1, psi_atom_indices + 1)
    phi_found_indices, phi_atom_indices = _atom_sequence(top, PHI_ATOMS)
    phi_info = _zip_dict(phi_found_indices + 1, phi_atom_indices + 1)
    dihedral_angle = {}
    for residue in top.residues:
        if residue.is_protein:
            dihedral_angle[residue.index+1] = {}
            if residue.index in phi_found_indices:
                dihedral_angle[residue.index+1]["phi"] = phi_info[residue.index+1]
            if residue.index in psi_found_indices:
                dihedral_angle[residue.index+1]["psi"] = psi_info[residue.index+1]
    return dihedral_angle


def get_dihedral_from_resid(file_path: str, selected_resid: List[int]) -> Dict:
    if len(selected_resid) == 0:
        return {}
    traj = md.load(file_path)
    top = traj.topology
    psi_found_indices, psi_atom_indices = _atom_sequence(top, PSI_ATOMS)
    psi_info = _zip_dict(psi_found_indices + 1, psi_atom_indices + 1)
    phi_found_indices, phi_atom_indices = _atom_sequence(top, PHI_ATOMS)
    phi_info = _zip_dict(phi_found_indices + 1, phi_atom_indices + 1)
    selected_dihedral_angle = {}
    residue_list = list(top.residues)
    for sid in selected_resid:
        residue = residue_list[sid-1]
        if residue.is_protein:
            selected_dihedral_angle[residue.index+1] = {}
            if residue.index in phi_found_indices:
                selected_dihedral_angle[residue.index+1]["phi"] = phi_info[residue.index+1]
            if residue.index in psi_found_indices:
                selected_dihedral_angle[residue.index+1]["psi"] = psi_info[residue.index+1]
    num_cv = len(selected_dihedral_angle.keys())
    logger.info(f"{num_cv} CVs have been created.")
    return selected_dihedral_angle

def get_distance_from_atomid(file_path: str, selected_atomid: List[int]) -> Dict:
    if len(selected_atomid) == 0:
        return {}
    top = md.load(file_path)
    selected_distance = {}
    for sid in selected_atomid:
        assert len(sid) == 2, "No valid distance list created."
        d_cv = md.compute_distances(top,atom_pairs=np.array([sid[0]-1,sid[1]-1]).reshape(-1,2),periodic=True)[0][0]
        selected_distance["%s %s"%(sid[0],sid[1])] = d_cv
    num_cv = len(selected_distance.keys())
    logger.info(f"{num_cv} CVs have been created.")
    return selected_distance

def slice_xtc_mdtraj(
        xtc: str,
        top: str,
        walker_idx: int,
        selected_idx: Sequence,
        output_format: str
    ):
    logger.info("slicing trajectories ...")
    traj = md.load_xtc(xtc, top=top)
    for sel in selected_idx:
        frame = traj[sel]
        frame.save_gro(output_format.format(walker=walker_idx,idx=sel))


def slice_xtc(
        xtc: str,
        top: str,
        walker_idx: int,
        selected_idx,
        output: str,
        style: str =  "gmx"
    ):
    if style == "gmx":
        slice_trjconv(
            xtc = xtc,
            top = top,
            selected_time = selected_idx,
            output = output
        )
    elif style == "mdtraj":
        slice_xtc_mdtraj(
            xtc = xtc,
            top = top,
            walker_idx=walker_idx,
            selected_idx = selected_idx,
            output_format=output
        )
    else:
        raise RuntimeError("Unknown Style for Slicing Trajectory.")