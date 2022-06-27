from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact
)

import json, shutil
from typing import Tuple, List, Optional, Dict
from pathlib import Path
from rid.constants import expore_task_pattern, md_mdp_name, plumed_input_name
from rid.common.gromacs import make_md_mdp_from_config
from rid.utils import write_txt


class PrepExplore(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "models": Optional[Artifact(List[Path])],
                "trust_level_0": float,
                "trust_level_1": float,
                "topology": Artifact(Path),
                "conformations": Artifact(List[Path]),
                "md_parameters": Dict,
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "msg": str,
                "task_list": Artifact(List[Path]),
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:
        assert op_in["topology"].exists()
        for idx, conf in enumerate(op_in["conformations"]):
            assert conf.exists()
            sub_task_path = Path(expore_task_pattern.format(idx))
            sub_task_path.mkdir(exist_ok=True, parents=True)
            shutil.copy(conf, sub_task_path.joinpath(conf.name))
            shutil.copy(op_in["topology"], sub_task_path.joinpath(op_in["topology"].name))
            sub_task_path.joinpath(md_mdp_name).write_text(
                make_md_mdp_from_config(op_in["md_parameters"])
            )
            sub_task_path.joinpath(md_mdp_name).write_text(
                make_md_mdp_from_config(op_in["md_parameters"])
            )


        out_msg = op_in["msg"]
        op_out = OPIO(
            {
                "msg": out_msg,
                "out_art": Path("bar.txt"),
            }
        )
        return op_out