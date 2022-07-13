from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Parameter
)

from typing import Tuple, List, Optional, Dict, Union
from pathlib import Path
from rid.select.cluster import Cluster
from rid.utils import save_txt, set_directory
from rid.constants import (
    culster_selection_data_name, 
    culster_selection_index_name, 
    cls_sel_precision, 
    cls_ndx_precision
)
import numpy as np


class PrepSelect(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "task_name": str,
                "plm_out": Artifact(Path),
                "cluster_threshold": float,
                "angular_mask": Optional[Union[np.ndarray, List]],
                "weights": Optional[Union[np.ndarray, List]],
                "numb_cluster_upper": Parameter(Optional[float], default=None),
                "numb_cluster_lower": Parameter(Optional[float], default=None),
                "max_selection": int,
                "if_make_threshold": Parameter(bool, default=False)
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "numb_cluster": int,
                "cluster_threshold": float,
                "culster_selection_index": Artifact(Path),
                "culster_selection_data": Artifact(Path)
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:
        # the first column of plm_out is time index
        data = np.loadtxt(op_in["plm_out"])[:,1:]
        cv_cluster = Cluster(
            data, op_in["cluster_threshold"], angular_mask=op_in["angular_mask"], 
            weights=op_in["weights"], max_selection=op_in["max_selection"])
        if op_in["if_make_threshold"]:
            assert (op_in["numb_cluster_lower"] is not None) and (op_in["numb_cluster_upper"] is not None), \
                "Please provide a number interval to make cluster thresholds."
            threshold = cv_cluster.make_threshold(op_in["numb_cluster_lower"], op_in["numb_cluster_upper"])
        else:
            threshold = op_in["cluster_threshold"]
        cls_sel_idx = cv_cluster.get_cluster_selection()
        selected_data = data[cls_sel_idx]
        numb_cluster = len(cls_sel_idx)

        task_path = Path(op_in["task_name"])
        task_path.mkdir(exist_ok=True, parents=True)
        with set_directory(task_path):
            save_txt(culster_selection_index_name, cls_sel_idx, fmt=cls_ndx_precision)
            save_txt(culster_selection_data_name, selected_data, fmt=cls_sel_precision)
        
        op_out = OPIO({
                "cluster_threshold": threshold,
                "numb_cluster": numb_cluster,
                "culster_selection_index": task_path.joinpath(culster_selection_index_name),
                "culster_selection_data": task_path.joinpath(culster_selection_data_name)
            })
        return op_out

    