from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Parameter
)

import json, shutil
from typing import Tuple, List, Optional, Dict
from pathlib import Path
from rid.utils import load_txt, save_txt
from rid.constants import sel_gro_name, cv_init_label, model_devi_name, model_devi_precision, sel_ndx_name
from rid.common.gromacs.command import get_grompp_cmd, get_mdrun_cmd
from rid.select.conf_select import select_from_devi
from rid.common.mol import slice_xtc
from rid.select.model_devi import make_std


class RunSelect(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "culster_selection_index": Artifact(Path),
                "culster_selection_data": Artifact(Path),
                "model_list": Artifact(List[Path]),
                "trust_lvl_1": float,
                "trust_lvl_2": float,
                "numb_cluster_threshold": float,
                "xtc_traj": Artifact(Path),
                "topology": Artifact(Path),
                "dt": Parameter(Optional[float], default=None),
                "slice_mode": Parameter(str, default="gmx")
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "selected_confs": Artifact(List[Path]),
                "selected_cv_init": Artifact(List[Path]),
                "number_of_cluster": int,
                "model_devi": Artifact(Path),
                "selected_indices": Artifact(Path)
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:
        cls_sel_idx = load_txt(op_in["culster_selection_index"], dtype=int)
        cls_sel_data = load_txt(op_in["culster_selection_data"], dtype=float)
        stds = make_std(cls_sel_data, models=op_in["model_list"])
        save_txt("cls_"+model_devi_name, stds, fmt=model_devi_precision)
        _selected_idx = select_from_devi(stds, op_in["trust_lvl_1"])
        sel_idx = cls_sel_idx[_selected_idx]
        save_txt(sel_ndx_name, sel_idx, fmt="%d")
        sel_data = cls_sel_data[_selected_idx]
        if op_in["slice_mode"] == "gmx":
            assert op_in["dt"] is not None, "Please provide time step to slice trajectory."
            for ii, sel in enumerate(sel_idx):
                time = sel * op_in["dt"]
                slice_xtc(xtc=op_in["xtc_traj"], top=op_in["topology"],
                        selected_idx=time, output=sel_gro_name.format(idx=sel), style="gmx")
        elif op_in["slice_mode"] == "mdtraj":
            slice_xtc(xtc=op_in["xtc_traj"], top=op_in["topology"],
                    selected_idx=sel_idx, output=sel_gro_name, style="mdtraj")
        else:
            raise RuntimeError("Unknown Style for Slicing Trajectory.")
        gro_list = []
        cv_init_list = []
        for ii, sel in enumerate(sel_idx):
            gro_list.append(Path(sel_gro_name.format(idx=sel)))
            save_txt(cv_init_label.format(idx=sel), sel_data[ii])
            cv_init_list.append(Path(cv_init_label.format(idx=sel)))
        numb_cluster = len(cls_sel_idx)
            
        op_out = OPIO(
            {
               "selected_confs": gro_list,
               "selected_cv_init": cv_init_list,
               "number_of_cluster": numb_cluster,
               "model_devi": Path(model_devi_name),
               "selected_indices": Path(sel_ndx_name)
            }
        )
        return op_out