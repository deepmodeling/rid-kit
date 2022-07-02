from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact
)

import json, shutil
from typing import Tuple, List, Optional, Dict
from pathlib import Path
from rid.constants import (
        label_task_pattern, 
        gmx_conf_name,
        gmx_top_name,
        gmx_mdp_name, 
        gmx_tpr_name,
        plumed_input_name,
        plumed_output_name,
        gmx_grompp_log,
        gmx_mdrun_log,
        trr_name,
        xtc_name
    )
from rid.utils import run_command, set_directory
from rid.common.gromacs.command import get_grompp_cmd, get_mdrun_cmd


class RunLabel(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "task_path", Artifact(Path)
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "msg": str,
                "task_files": Artifact(Path),
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:
        gmx_grompp_cmd = get_grompp_cmd(
            mdp = gmx_mdp_name,
            conf = gmx_conf_name,
            topology = gmx_top_name,
            output = gmx_tpr_name
        )
        gmx_run_cmd = get_mdrun_cmd(
            tpr=gmx_tpr_name,
            plumed=plumed_input_name
        )
        with set_directory(op_in["task_path"]):
            run_command(gmx_grompp_cmd)
            run_command(gmx_run_cmd)
        
        op_out = OPIO(
            {
                "plm_out": op_in["task_path"].joinpath(plumed_output_name),
                "gmx_grompp_log": op_in["task_path"].joinpath(gmx_grompp_log),
                "gmx_mdrun_log": op_in["task_path"].joinpath(gmx_mdrun_log)
            }
        )
        return op_out