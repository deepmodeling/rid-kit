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
        explore_task_pattern, 
        gmx_conf_name,
        gmx_top_name,
        gmx_mdp_name, 
        plumed_input_name,
        plumed_output_name
    )


class RunExplore(OP):

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
            plumed = plumed_input_name
        )
        gmx_run_cmd = get_mdrun_cmd(
            
        )
        
        op_out = OPIO(
            {
                "task_path": task_path,
            }
        )
        return op_out