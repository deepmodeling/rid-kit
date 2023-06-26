import dpdata
from typing import List, Dict, Sequence, Union

def slice_dump_dpdata(
    dump: str,
    walker_idx,
    selected_idx,
    output_format:str,
    type_map: List = [],
):
    if type_map != []:
        traj = dpdata.System(dump, fmt="lammps/dump", type_map=type_map)
    else:
        traj = dpdata.System(dump, fmt="lammps/dump")
    if output_format.endswith("lmp"):
        for sel in selected_idx:
            traj.to('lammps/lmp', output_format.format(walker=walker_idx,idx=sel), frame_idx=sel)
    elif output_format.endswith("gro"):
        for sel in selected_idx:
            traj.to('gromacs/gro', output_format.format(walker=walker_idx,idx=sel), frame_idx=sel)
    else:
        raise ValueError("Invalid output format, only support gmx and lmp")

def slice_dump(
    dump: str,
    walker_idx: int,
    selected_idx,
    output:str,
    type_map: List = [],
    style: str = "dpdata"
):
    if style == "dpdata":
        slice_dump_dpdata(
            dump = dump,
            walker_idx= walker_idx,
            selected_idx = selected_idx,
            type_map = type_map,
            output_format=output
        )
    else:
        raise RuntimeError("Unknown Style for Slicing Trajectory.")