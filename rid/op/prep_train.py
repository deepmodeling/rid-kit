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
from rid.task.builder import RestrainedMDTaskBuilder
from rid.utils import load_txt


class PrepExplore(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "topology": Artifact(Path),
                "conf": Artifact(Path),
                "gmx_config": Dict,
                "cv_config": Dict,
                "task_id": int,
                "kappas": List[float],
                "at": Artifact(List[Path()])
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
        if op_in["cv_config"]["mode"] == "torsion":
            cv_file = None,
            selected_resid = op_in["cv_config"]["selected_resid"]
        elif op_in["cv_config"]["mode"] == "custom":
            cv_file = op_in["cv_config"]["cv_file"],
            selected_resid = None
        at = load_txt(op_in["at"])
        gmx_task_builder = RestrainedMDTaskBuilder(
            conf = op_in["conf"],
            topology = op_in["topology"],
            gmx_config = op_in["gmx_config"],
            cv_file = cv_file,
            selected_resid = selected_resid,
            kappa = op_in["kappa"],
            at = at,
            plumed_output = plumed_output_name,
            cv_mode = op_in["cv_config"]["mode"]
        )
        gmx_task = gmx_task_builder.build()
        task_path = Path(explore_task_pattern.format(op_in["task_id"]))
        task_path.mkdir(exist_ok=True, parents=True)
        for fname, fconts in gmx_task.files.items():
            with open(task_path.joinpath(fname), fconts[1]) as ff:
                ff.write(fconts[0])
        op_out = OPIO(
            {
                "task_path": task_path
            }
        )
        return op_out