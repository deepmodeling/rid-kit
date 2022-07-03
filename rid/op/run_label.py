from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact
)
import os, sys
import logging
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
from rid.utils import run_command, set_directory, list_to_string
from rid.common.gromacs.command import get_grompp_cmd, get_mdrun_cmd

logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


class RunLabel(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "task_path": Artifact(Path),
                "md_config": Dict
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "plm_out": Artifact(Path),
                "gmx_grompp_log": Artifact(Path),
                "gmx_mdrun_log": Artifact(Path)
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
            output = gmx_tpr_name,
            max_warning=op_in["md_config"]["max_warning"]
        )
        gmx_run_cmd = get_mdrun_cmd(
            tpr=gmx_tpr_name,
            plumed=plumed_input_name,
            max_warning=op_in["md_config"]["max_warning"],
            nt=op_in["md_config"]["nt"],
            ntmpi=op_in["md_config"]["ntmpi"]
        )
        with set_directory(op_in["task_path"]):
            logger.info(list_to_string(gmx_grompp_cmd, " "))
            return_code, out, err = run_command(gmx_grompp_cmd)
            assert return_code == 0, err
            logger.info(err)

            logger.info(list_to_string(gmx_run_cmd, " "))
            return_code, out, err = run_command(gmx_run_cmd)
            assert return_code == 0, err
            logger.info(err)
        
        op_out = OPIO(
            {
                "plm_out": op_in["task_path"].joinpath(plumed_output_name),
                "gmx_grompp_log": op_in["task_path"].joinpath(gmx_grompp_log),
                "gmx_mdrun_log": op_in["task_path"].joinpath(gmx_mdrun_log)
            }
        )
        return op_out