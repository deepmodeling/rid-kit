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
from rid.utils import save_txt
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
                "plm_out": Artifact(Path),
                "cluster_threshold": float,
                "angular_mask": Optional[Union[np.ndarray, List]],
                "weights": Optional[Union[np.ndarray, List]],
                "numb_cluster_upper": float,
                "numb_cluster_lower": float,
                "max_selection": int,
                "if_make_threshold": Parameter(bool, default=False)
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "number_of_cluster": int,
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
            threshold = cv_cluster.make_threshold(op_in["numb_cluster_lower"], op_in["numb_cluster_upper"])
        else:
            threshold = op_in["cluster_threshold"]
        cls_sel_idx = cv_cluster.get_cluster_selection()
        save_txt(culster_selection_index_name, cls_sel_idx, fmt=cls_ndx_precision)
        selected_data = data[cls_sel_idx]
        save_txt(culster_selection_data_name, selected_data, fmt=cls_sel_precision)
        numb_cluster = len(cls_sel_idx)
        op_out = OPIO({
                "cluster_threshold": threshold,
                "number_of_cluster": numb_cluster,
                "culster_selection_index": Path(culster_selection_index_name),
                "culster_selection_data": Path(culster_selection_data_name)
            })
        return op_out

    