from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Parameter
)

from typing import List, Optional, Union
from pathlib import Path
from rid.select.cluster import Cluster
from rid.utils import save_txt, set_directory
from rid.constants import (
    cluster_selection_data_name, 
    cluster_selection_index_name,
    cluster_fig
)
import numpy as np


class PrepSelect(OP):

    """PrepSelect OP clusters CV outputs of each parallel walker from exploration steps and prepares representative 
    frames of each clusters for further selection steps.
    RiD-kit employs agglomerative clustering algorithm performed by Scikit-Learn python package. The distance matrix of CVs
    is pre-calculated, which is defined by Euclidean distance in CV space. For each cluster, one representive frame will 
    be randomly chosen from cluster members.
    For periodic collective variables, RiD-kit uses `angular_mask` to identify them and handle their periodic conditions 
    during distance calculation.
    In the first run of RiD iterations, PrepSelect will make a cluster threshold automatically from the initial guess of this value 
    and make cluter numbers of each parallel walker fall into the interval of `[numb_cluster_lower, numb_cluster_upper]`.
    """

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
                "cluster_fig": Artifact(Path, archive = None),
                "cluster_selection_index": Artifact(Path, archive = None),
                "cluster_selection_data": Artifact(Path, archive = None)
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

            - `task_name`: (`str`) Task names, used to make sub-directory for tasks.
            - `plm_out`: (`Artifact(Path)`) Outputs of CV values (`plumed.out` by default) from exploration steps.
            - `cluster_threshold`: (`float`) Cluster threshold of agglomerative clustering algorithm
            - `angular_mask`: (`array_like`) Mask for periodic collective variables. 1 represents periodic, 0 represents non-periodic.
            - `weights`: (`array_like`) Weights to cluster collective variables. see details in cluster parts.
            - `numb_cluster_upper`: (`Optional[float]`) Upper limit of cluster number to make cluster threshold.
            - `numb_cluster_lower`: (`Optional[float]`) Lower limit of cluster number to make cluster threshold.
            - `max_selection`: (`int`) Max selection number of clusters in Selection steps for each parallel walker.
                For each cluster, one representive frame will be randomly chosen from cluster members.
            - `if_make_threshold`: (`bool`) whether to make threshold to fit the cluster number interval. Usually `True` in the 1st 
                iteration and `False` in the further iterations. 

        Returns
        -------
            Output dict with components:
        
            - `numb_cluster`: (`int`) Number of clusters.
            - `cluster_threshold`: (`float`) Cluster threshold of agglomerative clustering algorithm. 
            - `cluster_selection_index`: (`Artifact(Path)`) Indice of chosen representive frames of clusters in trajectories.
            - `cluster_selection_data`: (`Artifact(Path)`) Collective variable values of chosen representive frames of clusters.
        """
        task_path = Path(op_in["task_name"])
        task_path.mkdir(exist_ok=True, parents=True)
        # the first column of plm_out is time index, the second columnn is biased potential
        data = np.loadtxt(op_in["plm_out"])[:,2:]
        cv_cluster = Cluster(
            data, op_in["cluster_threshold"], op_in["task_name"], angular_mask=op_in["angular_mask"], 
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

        with set_directory(task_path):
            np.save(cluster_selection_index_name, cls_sel_idx)
            np.save(cluster_selection_data_name, selected_data)
        
        op_out = OPIO({
                "cluster_threshold": threshold,
                "numb_cluster": numb_cluster,
                "cluster_fig": task_path.joinpath(cluster_fig),
                "cluster_selection_index": task_path.joinpath(cluster_selection_index_name),
                "cluster_selection_data": task_path.joinpath(cluster_selection_data_name)
            })
        return op_out

    