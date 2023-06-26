from typing import List, Optional, Dict
from pathlib import Path
import numpy as np
from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Parameter,
    BigParameter
)
from rid.utils import save_txt, set_directory
from rid.constants import sel_gro_name, cv_init_label, model_devi_name, model_devi_precision, sel_ndx_name, new_model_devi_fig
from rid.select.conf_select import select_from_devi
from rid.common.mol import get_distance_from_atomid
from rid.select.model_devi import make_std
import os
from matplotlib import pyplot as plt

class RunModelDevi(OP):

    """
    `RunModelDevi` calculates model deviations for each chosen representive cluster frames from `PrepSelect` and select 
    ones with high uncertainty from them.
    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "models": Artifact(List[Path]),
                "plm_out": Artifact(Path),
                "selected_indices": Artifact(Path, optional=True),
                "trust_lvl_1": float,
                "task_name": str,
                "block_tag": str
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "model_devi": Artifact(Path, optional=True, archive = None),
                "model_devi_fig": Artifact(Path, optional=True, archive = None)
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

            - `models`: Artifact(List[Path], optional=True),
            - `plm_out`: Artifact(Path)
          
        Returns
        -------
            Output dict with components:
        
            - `model_devi`: (`Artifact(Path)`) Model deviation values of selected conformation files (`selected_confs`) .
        """
        task_path = Path(op_in["task_name"])
        task_path.mkdir(exist_ok=True, parents=True)

        with set_directory(task_path):
            # if "units" in cv_config:
            #     length_units = cv_config["units"]
            # else:
            #     length_units = None
            # if length_units == None or length_units == "nm":
            #     cv_data = np.loadtxt(op_in["plm_out"])[:,2:]
            # elif length_units == "A" or length_units == "Angstrom":
            #     cv_data = 10*np.loadtxt(op_in["plm_out"])[:,2:]
            # else:
            #     raise ValueError("Valid length units must be nm or A")
            trust_lvl = op_in["trust_lvl_1"]
            cv_data = np.loadtxt(op_in["plm_out"])[:,2:]
            stds = make_std(cv_data, models=op_in["models"])
            save_txt("cls_"+model_devi_name, stds, fmt=model_devi_precision)
            
            _selected_idx = np.load(op_in["selected_indices"])
     
            nframes = cv_data.shape[0]
            groups = ["other"]*nframes
            xlist = np.linspace(0,nframes-1,nframes,dtype=int)
            ylist = [int(j) for j in _selected_idx]                
            for y_ in ylist:
                groups[y_] = "train"
            cdict = {"other": 'black',"train": 'red'}
            fig, ax = plt.subplots()
            groups = np.array(groups)
            
            for g in ["other","train"]:
                ix = np.where(groups == g)
                ax.scatter(xlist[ix], stds[ix], c = cdict[g], label = g, s = 20)
            plt.axhline(y=trust_lvl, color='b', linestyle='-',label="trust level 1")
            plt.xlabel("configuration list")
            plt.ylabel("force std (KJ/(mol*length))")
            plt.legend()
            plt.title("prediction error")
            fig.set_size_inches(8, 6, forward=True)
            plt.savefig(new_model_devi_fig)
                
        model_devi_fig_path = None
        if os.path.exists(task_path.joinpath(new_model_devi_fig)):
            model_devi_fig_path = task_path.joinpath(new_model_devi_fig)
            
        op_out = OPIO(
            {
               "model_devi": task_path.joinpath("cls_"+model_devi_name),
               "model_devi_fig": model_devi_fig_path
            }
        )
        return op_out