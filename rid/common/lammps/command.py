import dpdata

def final_dump(
    dump: str,
    selected_idx,
    output_format:str
):
    traj = dpdata.System(dump, fmt="lammps/dump")
    traj.to('lammps/lmp', output_format, frame_idx=selected_idx)