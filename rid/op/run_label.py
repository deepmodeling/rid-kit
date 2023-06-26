import os, sys
import logging
from typing import Dict, List
from pathlib import Path
from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Parameter,
    BigParameter
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
        gmx_xtc_name,
        gmx_trr_name,
        gmx_coord_name,
        gmx_force_name,
        mf_fig
    )

from rid.utils import run_command, set_directory, list_to_string
from rid.common.sampler.command import get_grompp_cmd, get_mdrun_cmd
from rid.common.gromacs.trjconv import generate_coords, generate_forces
from rid.tools.estimator import CalcMF
import numpy as np

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
    where harmonnic restraints are exerted on collective variables.
    `RunLabel` is able to run in a standard Gromacs-PLUMED2 or Lammps-PLUMED2 env.
    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "task_path": Artifact(Path),
                "forcefield": Artifact(Path, optional=True),
                "label_config": BigParameter(Dict),
                "cv_config": BigParameter(Dict),
                "at": Artifact(Path, optional=True),
                "task_name": BigParameter(str),
                "tail": Parameter(float, default=0.9),
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
                "plm_out": Artifact(Path, archive = None),
                "trajectory": Artifact(Path, archive = None),
                "cv_forces": Artifact(Path, archive = None),
                "mf_info": Artifact(Path,optional=True, archive = None),
                "mf_fig": Artifact(Path, archive = None),
                "md_log": Artifact(Path, archive = None)
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

            - `task_path`: (`Artifact(Path)`) A directory path containing files for Gromacs/Lammps MD simulations.
            - `label_config`: (`Dict`) Configuration of Gromacs/Lammps simulations in label steps.
            - `cv_config`: (`Dict`) Configuration of CV in MD simulations.
          
        Returns
        -------
            Output dict with components:
        
            - `forces`: (`Artifact(Path)`) mean forces value of restrained/constrained estimator.
            - `mf_info`: (`Artifact(Path)`) mean force value and std information of restrained/constrained estimator.
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
                if not os.path.islink(op_in["forcefield"].name):
                    os.symlink(op_in["forcefield"], op_in["forcefield"].name)
            if op_in["dp_files"] is not None:
                for file in op_in["dp_files"]:
                    if not os.path.islink(file.name):
                        os.symlink(file, file.name)
            if op_in["index_file"] is not None:
                if not os.path.islink(op_in["index_file"].name):
                    os.symlink(op_in["index_file"], op_in["index_file"].name)
            if op_in["cv_file"] is not None:
                for file in op_in["cv_file"]:
                    if file.name != "colvar":
                        if not os.path.islink(file.name):
                            os.symlink(file, file.name)
            if op_in["inputfile"] is not None:
                for file in op_in["inputfile"]:
                    if file.name == inputfile_name:
                        if not os.path.islink(file.name):
                            os.symlink(file, file.name)
            
            if op_in["dp_files"] is not None:
                for dp_file in op_in["dp_files"]:
                    if (dp_file.name).endswith("json"):
                        os.environ["GMX_DEEPMD_INPUT_JSON"] = dp_file.name
                
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
            
            if op_in["label_config"]["method"] == "constrained":
                select_system = "System"
                if "system" in op_in["label_config"]:
                    select_system = op_in["label_config"]["system"]
                generate_coords(output_group = select_system, trr = gmx_trr_name, top = op_in["task_path"].joinpath(gmx_conf_name), index = op_in["index_file"], out_coord=gmx_coord_name)
                generate_forces(output_group = select_system, trr = gmx_trr_name, top = op_in["task_path"].joinpath(gmx_conf_name), index = op_in["index_file"], out_force=gmx_force_name)

        frame_coords = None
        frame_forces = None
        traj_path = None
        if op_in["label_config"]["type"] == "gmx":
            mdrun_log = gmx_mdrun_log
            if "traj_output" in op_in["label_config"] and op_in["label_config"]["traj_output"] == "True":
                if os.path.exists(op_in["task_path"].joinpath(gmx_xtc_name)):
                    traj_path = op_in["task_path"].joinpath(gmx_xtc_name)
            if op_in["label_config"]["method"] == "constrained":
                frame_coords = op_in["task_path"].joinpath(gmx_coord_name)
                frame_forces = op_in["task_path"].joinpath(gmx_force_name)
        elif op_in["label_config"]["type"] == "lmp":
            mdrun_log = lmp_mdrun_log
            if "traj_output" in op_in["label_config"] and op_in["label_config"]["traj_output"] == "True":
                if os.path.exists(op_in["task_path"].joinpath("out.dump")):
                    traj_path = op_in["task_path"].joinpath("out.dump")
        
        conf = None
        topology = None
        plm_out = None
        
        if os.path.exists(op_in["task_path"].joinpath(gmx_conf_name)):
            conf = op_in["task_path"].joinpath(gmx_conf_name)
        if os.path.exists(op_in["task_path"].joinpath(gmx_top_name)):
            topology = op_in["task_path"].joinpath(gmx_top_name)
        if os.path.exists(op_in["task_path"].joinpath(plumed_output_name)):
            plm_out = op_in["task_path"].joinpath(plumed_output_name)
            
        mf_out = CalcMF(
            conf = conf,
            task_name = op_in["task_name"],
            plm_out = plm_out,
            cv_config = op_in["cv_config"],
            label_config = op_in["label_config"],
            tail = op_in["tail"],
            topology = topology,
            frame_coords = frame_coords,
            frame_forces = frame_forces,
            at = op_in["at"]
        )
            
        op_out = OPIO(
            {
                "plm_out": plm_out,
                "trajectory": traj_path,
                "cv_forces": mf_out["cv_forces"],
                "mf_info":  mf_out["mf_info"],
                "mf_fig": Path(op_in["task_name"]).joinpath(mf_fig),
                "md_log": op_in["task_path"].joinpath(mdrun_log)
            }
        )
        return op_out