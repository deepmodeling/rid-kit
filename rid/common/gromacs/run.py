from typing import Optional, Union, List
from rid.common.gromacs.gmx_constant import gmx_prep_cmd

def get_grompp_cmd(
        mdp: Optional[str] = None,
        conf: Optional[str] = None,
        topology: Optional[str] = None,
        ref: Optional[str] = None,
        index: Optional[str] = None,
        plumed: Optional[str] = None,
        max_warning: Optional[int] = None,
        output: Optional[str] = None,
        extra_parameters: Optional[Union[List, str]] = None
    ):
    cmd_list = [gmx_prep_cmd]
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
    if plumed is not None:
        cmd_list += ["-p", plumed]
    if max_warning is not None:
        cmd_list += ["-maxwarn", conf]
    if output is not None:
        cmd_list += ["-o", output]
    if extra_parameters is not None:
        if isinstance(extra_parameters, List):
            cmd_list += extra_parameters
        elif isinstance(extra_parameters, str):
            cmd_list.append(extra_parameters)
        else:
            raise RuntimeError("Unknown extra parameters format.")
    return " ".join(cmd_list)

        



def mdrun_cmd():

    return 