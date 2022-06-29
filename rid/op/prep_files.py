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
        md_mdp_name, 
        plumed_input_name
    )
from rid.common.gromacs import make_md_mdp_from_config
from rid.common.plumed import make_deepfe_plumed, check_deepfe_input
from rid.utils import write_txt, load_txt, load_pkl


class PrepExplore(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "task": Artifact(Path),
                "task_id": int,
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "task_path": Artifact(Path),
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:
        task = load_pkl(op_in["task"])
        for fname, fconts in  task.files.items():
            with open(fname, fconts[1]) as ff:
                ff.write(fconts[0])
        op_out = OPIO(
            {
                "task_path": explore_task_pattern.format(op_in["task_id"]),
            }
        )
        return op_out