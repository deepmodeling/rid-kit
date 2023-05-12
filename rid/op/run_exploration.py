import os, sys
import logging
from typing import List, Dict, Union
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
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
        gmx_conf_out,
        lmp_conf_out,
        bias_fig,
        model_devi_fig,
        dp_model_devi_fig,
        dp_sel_ndx,
        sel_lmp_name,
        lmp_input_name
    )
from rid.utils import run_command, set_directory, list_to_string
from rid.common.lammps.command import final_dump
from rid.common.sampler.command import get_grompp_cmd, get_mdrun_cmd
from rid.common.mol_dpdata import slice_dump
from rid.select.model_devi import make_std
from rid.mcmc.walker import s,z


logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


class RunExplore(OP):

    """Run (biased or brute-force) MD simulations with files provided by `PrepExplore` OP.
    RiD-kit emploies Gromacs/Lammps as MD engine with PLUMED2 plugin (with or without MPI-implement). Make sure work environment 
    has properly installed Gromacs/Lammps with patched PLUMED2. Also, PLUMED2 need an additional `DeePFE.cpp` patch to use bias potential of RiD.
    See `install` in the home page for details.
    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "task_path": Artifact(Path),
                "forcefield": Artifact(Path, optional=True),
                "exploration_config": BigParameter(Dict),
                "models": Artifact(List[Path], optional=True),
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
                "bias_fig": Artifact(Path, optional=True, archive = None),
                "model_devi_fig": Artifact(Path,optional=True, archive = None),
                "dp_model_devi_fig": Artifact(Path,optional=True, archive = None),
                "dp_model_devi": Artifact(Path,optional=True, archive = None),
                "dp_selected_indices": Artifact(Path,optional=True, archive = None),
                "dp_selected_confs": Artifact(List[Path], archive = None),
                "projected_fig": Artifact(Path,optional=True, archive = None),
                "md_log": Artifact(Path, archive = None),
                "trajectory": Artifact(Path, archive = None),
                "conf_out": Artifact(Path, archive = None)
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
            - `exploration_config`: (`Dict`) Configuration of Gromacs/Lammps simulations in exploration steps.
            - `models`: (`Artifact(List[Path])`) Optional. Neural network model files (`.pb`) used to bias the simulation. 
                Run brute force MD simulations if not provided.
          
        Returns
        -------
            Output dict with components:
        
            - `plm_out`: (`Artifact(Path)`) Outputs of CV values (`plumed.out` by default) from exploration steps.
            - `md_log`: (`Artifact(Path)`) Log files of Gromacs/lammps `mdrun` commands.
            - `trajectory`: (`Artifact(Path)`) Trajectory files (`.xtc`). The output frequency is defined in `exploration_config`.
            - `conf_out`: (`Artifact(Path)`) Final frames of conformations in simulations.
        """
        if op_in["index_file"] is None:
            index = None
        else:
            index = op_in["index_file"].name
        
        if "max_warning" in op_in["exploration_config"]:
            max_warning=op_in["exploration_config"]["max_warning"]
        else:
            max_warning = None
        
        if "nt" in op_in["exploration_config"]:
            nt = op_in["exploration_config"]["nt"]
        else:
            nt = None
        
        if "ntmpi" in op_in["exploration_config"]:
            ntmpi = op_in["exploration_config"]["ntmpi"]
        else:
            ntmpi = None
        
        if "inputfile" in op_in["exploration_config"]:
            inputfile_name = op_in["exploration_config"]["inputfile"]
        else:
            inputfile_name = None
        
        grompp_cmd = get_grompp_cmd(
            sampler_type=op_in["exploration_config"]["type"],
            mdp = gmx_mdp_name,
            conf = gmx_conf_name,
            topology = gmx_top_name,
            output = gmx_tpr_name,
            index = index,
            max_warning=max_warning
        )
        run_cmd = get_mdrun_cmd(
            sampler_type=op_in["exploration_config"]["type"],
            tpr=gmx_tpr_name,
            plumed=plumed_input_name,
            nt=nt,
            ntmpi=ntmpi,
            inputfile=inputfile_name
        )

        with set_directory(op_in["task_path"]):
            if op_in["forcefield"] is not None:
                if not os.path.islink(op_in["forcefield"].name):
                    os.symlink(op_in["forcefield"], op_in["forcefield"].name)
            if op_in["models"] is not None:
                for model in op_in["models"]:
                    if not os.path.islink(model.name):
                        os.symlink(model, model.name)
            if op_in["dp_files"] is not None:
                for file in op_in["dp_files"]:
                    if not os.path.islink(file.name):
                        os.symlink(file, file.name)

            if op_in["cv_file"] is not None:
                for file in op_in["cv_file"]:
                    if file.name != "colvar":
                        if not os.path.islink(file.name):
                            os.symlink(file, file.name)
                            
            if op_in["exploration_config"]["type"] == "lmp":
                if op_in["inputfile"] is not None:
                    for file in op_in["inputfile"]:
                        if file.name == inputfile_name:
                            if not os.path.islink(file.name):
                                os.symlink(file, file.name)
            
            if op_in["exploration_config"]["type"] == "gmx":
                if op_in["index_file"] is not None:
                    if not os.path.islink(op_in["index_file"].name):
                        os.symlink(op_in["index_file"], op_in["index_file"].name)
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
            if op_in["exploration_config"]["type"] == "lmp":
                final_dump(dump="out.dump",selected_idx=-1,output_format=lmp_conf_out)
            
            bias_fig_file = None
            model_devi_fig_file = None
            
            dp_model_devi_fig_file = None
            dp_model_devi_file = None
            dp_sel_ndx_file = None
            dp_conf_list = []
            if "dp_model_devi_out" in op_in["exploration_config"]:
                # plot DP potential model deviation
                dp_model_devi_file_name = op_in["exploration_config"]["dp_model_devi_out"]
                data = np.loadtxt(dp_model_devi_file_name)
                f_list = data[:,4]
                cv_data = np.loadtxt("plm.out")[:,2:]
                assert (len(cv_data) == len(f_list))
                trust_lvl_hi = op_in["exploration_config"]["dp_model_devi_f_hi"]
                trust_lvl_lo = op_in["exploration_config"]["dp_model_devi_f_lo"]
                max_selection = op_in["exploration_config"]["dp_model_max_selection"]
                high_devi = np.where((f_list>trust_lvl_lo) & (f_list < trust_lvl_hi))
                
                if len(high_devi[0]) > max_selection:
                    high_devi = np.random.choice(high_devi[0],max_selection,replace=False)
                else:
                    high_devi = high_devi[0]
                np.savetxt(dp_sel_ndx,high_devi,fmt="%d")
                dp_sel_ndx_file = op_in["task_path"].joinpath(dp_sel_ndx)
                walker_idx = int(op_in["task_path"].name)
                slice_dump(dump="out.dump",walker_idx = walker_idx,selected_idx=high_devi, output=sel_lmp_name,style="dpdata")
                for ii, sel in enumerate(high_devi):
                    dp_conf_list.append(op_in["task_path"].joinpath(sel_lmp_name.format(walker=walker_idx,idx=sel)))
                
                x_list = np.array([i for i in range(len(f_list))])
                nframes = len(x_list)
                groups = ["other"]*nframes
                y_list = [int(j) for j in high_devi]
                for j in y_list:
                    groups[j] = "choice"
                cdict = {"other": 'black', "choice": 'red'}
                groups = np.array(groups)
                fig, ax = plt.subplots()
                for g in ["other","choice"]:
                    ix = np.where(groups == g)
                    ax.scatter(x_list[ix], f_list[ix], c = cdict[g], label = g, s = 20)
                plt.title("max model devi of DP model")
                plt.legend()
                plt.xlabel("simu frames")
                plt.ylabel("max model devi")
                plt.savefig(dp_model_devi_fig)
                dp_model_devi_file = op_in["task_path"].joinpath(dp_model_devi_file_name)
                dp_model_devi_fig_file = op_in["task_path"].joinpath(dp_model_devi_fig)
            
            projected_fig_file = None
            if "proj_info" in op_in["exploration_config"]:
                proj_info = op_in["exploration_config"]["proj_info"]
                # plot projected cv fig
                proj_mode = proj_info["proj_mode"]
                if proj_mode == "path":
                    proj_cv_index = proj_info["proj_cv_index"]
                    proj_lm = proj_info["path_lm"]
                    path_list = proj_info["path_list"]
                    cv_data = np.loadtxt("plm.out")[:,2:]
                    proj_cv_data = cv_data[:,proj_cv_index]                                      
                    slist = s(proj_cv_data, proj_lm, path_list)
                    zlist = z(proj_cv_data, proj_lm, path_list)
                    point_numbers = list(range(len(slist)))
                    plt.figure(figsize=(10, 8), dpi=100)
                    plt.scatter(slist,zlist, c = point_numbers)
                    plt.xlabel("s")
                    plt.ylabel("z")
                    plt.title("Trajectory projection")
                    plt.colorbar(label="steps")
                    plt.savefig("path.png")
                    projected_fig_file = op_in["task_path"].joinpath("path.png")
            
            if op_in["models"] is not None:
                # plot bias during simulation
                bias = np.loadtxt(plumed_output_name)[:,1]
                dt = op_in["exploration_config"]["dt"]*op_in["exploration_config"]["output_freq"]
                xlist = [i*dt for i in range(len(bias))]
                plt.figure(figsize=(10, 8), dpi=100)
                plt.scatter(xlist,np.log(bias+0.001))
                plt.xlabel("simulation time (ps)")
                plt.ylabel("log of bias (kj/(mol))")
                plt.title("log bias for the simulation")
                plt.savefig(op_in["task_path"].joinpath(bias_fig))
                # plot model deviation during simulation
                cv_list = np.loadtxt(plumed_output_name)[:,2:]
                stds = make_std(cv_list, op_in["models"])
                plt.figure(figsize=(10, 8), dpi=100)
                plt.scatter(xlist,stds)
                plt.xlabel("simulation time (ps)")
                plt.ylabel("force std (KJ/(mol*nm))")
                plt.title("model deviation during simulation")
                plt.savefig(op_in["task_path"].joinpath(model_devi_fig))
                
                bias_fig_file = op_in["task_path"].joinpath(bias_fig)
                model_devi_fig_file = op_in["task_path"].joinpath(model_devi_fig)
                
                        
        if op_in["exploration_config"]["type"] == "gmx":
            mdrun_log = gmx_mdrun_log
            traj_name = gmx_xtc_name
            conf_out = gmx_conf_out
        elif op_in["exploration_config"]["type"] == "lmp":
            mdrun_log = lmp_mdrun_log
            traj_name = "out.dump"
            conf_out = lmp_conf_out
        else:
            raise RuntimeError("Invalid sampler type")
            
        op_out = OPIO(
            {
                "plm_out": op_in["task_path"].joinpath(plumed_output_name),
                "bias_fig": bias_fig_file,
                "model_devi_fig": model_devi_fig_file,
                "dp_model_devi_fig": dp_model_devi_fig_file,
                "dp_model_devi": dp_model_devi_file,
                "dp_selected_indices": dp_sel_ndx_file,
                "dp_selected_confs": dp_conf_list,
                "projected_fig": projected_fig_file,
                "md_log": op_in["task_path"].joinpath(mdrun_log),
                "trajectory": op_in["task_path"].joinpath(traj_name),
                "conf_out": op_in["task_path"].joinpath(conf_out),
            }
        )
        return op_out