from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Parameter
)
import os, sys
import logging
import json, shutil
from typing import Tuple, List, Optional, Dict
from pathlib import Path
from rid.constants import (
        explore_task_pattern, 
        gmx_conf_name,
        gmx_top_name,
        gmx_mdp_name, 
        gmx_tpr_name,
        plumed_input_name,
        plumed_output_name,
        gmx_grompp_log,
        gmx_mdrun_log,
        xtc_name,
        gmx_conf_out
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


class RunExplore(OP):

    """Run (biased or brute-force) MD simulations with files provided by `PrepExplore` OP.
    RiD-kit emploies Gromacs as MD engine with PLUMED2 plugin (with or without MPI-implement). Make sure work environment 
    has properly installed Gromacs with patched PLUMED2. Also, PLUMED2 need an additional `DeePFE.cpp` patch to use bias potential of RiD.
    See `install` in the home page for details.
    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "task_path": Artifact(Path),
                "gmx_config": Dict,
                "models": Artifact(List[Path], optional=True)
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "plm_out": Artifact(Path),
                "md_log": Artifact(Path),
                "trajectory": Artifact(Path),
                "conf_out": Artifact(Path)
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:

        r"""Execute the OP.
        Parameters
        ----------
        op_in : dict
            Input dict with components:

            - `task_path`: (`Artifact(Path)`) A directory path containing files for Gromacs MD simulations.
            - `gmx_config`: (`Dict`) Configuration of Gromacs simulations in exploration steps.
            - `models`: (`Artifact(List[Path])`) Optional. Neural network model files (`.pb`) used to bias the simulation. 
                Run brute force MD simulations if not provided.
          
        Returns
        -------
            Output dict with components:
        
            - `plm_out`: (`Artifact(Path)`) Outputs of CV values (`plumed.out` by default) from exploration steps.
            - `md_log`: (`Artifact(Path)`) Log files of Gromacs `mdrun` commands.
            - `trajectory`: (`Artifact(Path)`) Trajectory files (`.xtc`). The output frequency is defined in `gmx_config`.
            - `conf_out`: (`Artifact(Path)`) Final frames of conformations in simulations.
        """

        gmx_grompp_cmd = get_grompp_cmd(
            mdp = gmx_mdp_name,
            conf = gmx_conf_name,
            topology = gmx_top_name,
            output = gmx_tpr_name,
            max_warning=op_in["gmx_config"]["max_warning"]
        )
        gmx_run_cmd = get_mdrun_cmd(
            tpr=gmx_tpr_name,
            plumed=plumed_input_name,
            nt=op_in["gmx_config"]["nt"],
            ntmpi=op_in["gmx_config"]["ntmpi"]
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
                "md_log": op_in["task_path"].joinpath(gmx_mdrun_log),
                "trajectory": op_in["task_path"].joinpath(xtc_name),
                "conf_out": op_in["task_path"].joinpath(gmx_conf_out),
            }
        )
        return op_out