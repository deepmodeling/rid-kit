import os, sys
import logging
from typing import Dict, List
from pathlib import Path
from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact
)
from rid.constants import (
        gmx_conf_name,
        gmx_top_name,
        gmx_mdp_name, 
        gmx_tpr_name,
        plumed_input_name,
        plumed_output_name,
        gmx_mdrun_log,
        lmp_mdrun_log,
        lmp_input_name
    )
from rid.utils import run_command, set_directory, list_to_string
from rid.common.sampler.command import get_grompp_cmd, get_mdrun_cmd

logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


class RunLabel(OP):

    """
    In `RunLabel`, labeling processes are achieved by restrained MD simulations 
    wehre harmonnic restraints are exerted on collective variables.
    `RunLabel` is able to run in a standard Gromacs-PLUMED2 env.
    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "task_path": Artifact(Path),
                "forcefield": Artifact(Path, optional=True),
                "label_config": Dict,
                "index_file": Artifact(Path, optional=True),
                "dp_files": Artifact(List[Path], optional=True),
                "cv_file": Artifact(List[Path], optional=True),
                "inputfile": Artifact(List[Path], optional=True)
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "plm_out": Artifact(Path),
                "md_log": Artifact(Path)
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
            - `label_config`: (`Dict`) Configuration of Gromacs simulations in label steps.
          
        Returns
        -------
            Output dict with components:
        
            - `plm_out`: (`Artifact(Path)`) Outputs of CV values (`plumed.out` by default) from label steps.
            - `md_log`: (`Artifact(Path)`) Log files of Gromacs `mdrun` commands.
        """
        
        if op_in["index_file"] is None:
            index = None
        else:
            index = op_in["index_file"].name
            
        if "inputfile" in op_in["label_config"]:
            inputfile_name = op_in["label_config"]["inputfile"]
        else:
            inputfile_name = None
            
        grompp_cmd = get_grompp_cmd(
            sampler_type=op_in["label_config"]["type"],
            mdp = gmx_mdp_name,
            conf = gmx_conf_name,
            topology = gmx_top_name,
            output = gmx_tpr_name,
            index = index,
            max_warning=op_in["label_config"]["max_warning"]
        )
        run_cmd = get_mdrun_cmd(
            sampler_type=op_in["label_config"]["type"],
            tpr=gmx_tpr_name,
            plumed=plumed_input_name,
            nt=op_in["label_config"]["nt"],
            ntmpi=op_in["label_config"]["ntmpi"],
            inputfile=inputfile_name
        )

        with set_directory(op_in["task_path"]):
            if op_in["forcefield"] is not None:
                os.symlink(op_in["forcefield"], op_in["forcefield"].name)
            if op_in["dp_files"] is not None:
                for file in op_in["dp_files"]:
                    os.symlink(file, file.name)
            if op_in["index_file"] is not None:
                os.symlink(op_in["index_file"], op_in["index_file"].name)
            if op_in["cv_file"] is not None:
                for file in op_in["cv_file"]:
                    if file.name != "colvar":
                        os.symlink(file, file.name)
            if op_in["inputfile"] is not None:
                for file in op_in["inputfile"]:
                    if file.name == inputfile_name:
                        os.symlink(file, file.name)
            
            if op_in["dp_files"] is not None:
                os.environ["GMX_DEEPMD_INPUT_JSON"] = "./input.json"
                
            if grompp_cmd is not None:
                logger.info(list_to_string(grompp_cmd, " "))
                return_code, out, err = run_command(grompp_cmd)
                assert return_code == 0, err
                logger.info(err)
            if run_cmd is not None:
                logger.info(list_to_string(run_cmd, " "))
                return_code, out, err = run_command(run_cmd)
                assert return_code == 0, err
                logger.info(err)
    
        if op_in["label_config"]["type"] == "gmx":
            mdrun_log = gmx_mdrun_log
        elif op_in["label_config"]["type"] == "lmp":
            mdrun_log = lmp_mdrun_log
            
        op_out = OPIO(
            {
                "plm_out": op_in["task_path"].joinpath(plumed_output_name),
                "md_log": op_in["task_path"].joinpath(mdrun_log)
            }
        )
        return op_out