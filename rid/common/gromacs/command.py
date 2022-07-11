from typing import Optional, Union, List
from rid.common.gromacs.gmx_constant import gmx_prep_cmd, gmx_run_cmd

def get_grompp_cmd(
        mdp: Optional[str] = None,
        conf: Optional[str] = None,
        topology: Optional[str] = None,
        ref: Optional[str] = None,
        index: Optional[str] = None,
        max_warning: Optional[int] = None,
        output: Optional[str] = None,
        extra_parameters: Optional[Union[List, str]] = None
    ):
    cmd_list = gmx_prep_cmd.split()
    if mdp is not None:
        cmd_list += ["-f", mdp]
    if conf is not None:
        cmd_list += ["-c", conf]
    if topology is not None:
        cmd_list += ["-p", topology]
    if ref is not None:
        cmd_list += ["-r", ref]
    if index is not None:
        cmd_list += ["-n", index]
    if max_warning is not None:
        cmd_list += ["-maxwarn", f"{max_warning:d}"]
    if output is not None:
        cmd_list += ["-o", output]
    if extra_parameters is not None:
        if isinstance(extra_parameters, List):
            cmd_list += extra_parameters
        elif isinstance(extra_parameters, str):
            cmd_list.append(extra_parameters)
        else:
            raise RuntimeError("Unknown extra parameters format.")
    return cmd_list


def get_mdrun_cmd(
        tpr: Optional[str],
        plumed: Optional[str] = None,
        cpi: Optional[str] = None,
        ntmpi: Optional[int] = None,
        nt: Optional[int] = None,
        extra_parameters: Optional[Union[List, str]] = None
    ):
    run_cmd = gmx_run_cmd.split()
    if tpr is not None:
        run_cmd += ["-s", tpr]
    if cpi is not None:
        run_cmd += ["-cpi", cpi]
    if plumed is not None:
        run_cmd += ["-plumed", plumed]
    if ntmpi is not None:
        run_cmd += ["-ntmpi", f"{ntmpi:d}"]
    if nt is not None:
        run_cmd += ["-nt", f"{nt:d}"]
    if extra_parameters is not None:
        if isinstance(extra_parameters, List):
            run_cmd += extra_parameters
        elif isinstance(extra_parameters, str):
            run_cmd.append(extra_parameters)
        else:
            raise RuntimeError("Unknown extra parameters format.")
    return run_cmd