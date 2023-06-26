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
import numpy as np
from rid.utils import set_directory
from matplotlib import pyplot as plt
from rid.constants import mf_std_fig
import os


class LabelStats(OP):

    r"""LabelStats OP calculate the std of all the labeling steps, remove the ones with high std.
    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "cv_forces": Artifact(List[Path]),
                "mf_info": Artifact(List[Path]), 
                "std_threshold": float  
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "mf_std_fig":  Artifact(Path, archive=None),
                "cv_forces": Artifact(List[Path], archive=None)
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
        
            - "cv_forces": (`Artifact(List[Path])`) 
            - "mf_info": (`Artifact(List[Path])`) 
            - "std_threshold": (`Float`)

        Returns
        -------
            Output dict with components:
        
            - "mf_std_fig":  (`Artifact(Path)`)
            - "cv_forces":(`Artifact(List[Path])`) 
            
        """
        assert len(op_in["mf_info"]) == len(op_in["cv_forces"])
        cv_forces_list = [_ for _ in op_in["cv_forces"] if _ is not None]
        cv_forces_list = np.array(cv_forces_list)
        cv_force_file = cv_forces_list[0]
        cv_force = np.loadtxt(cv_force_file)
        cv_dim = int(cv_force.shape[0]//2)
        # extract the std from the mf_info
        mf_all_std_list = []
        for mf_info in op_in["mf_info"]:
            if mf_info and os.path.exists(mf_info):
                with open(mf_info) as f:
                    while True:
                        line = f.readline()
                        line_list = line.strip().split(" ")[:3]
                        if line_list == ['mean', 'force', 'std']:
                            break
                    mf_std_line = line.strip().split(" ")[7:]
                    mf_std_list = [float(i) for i in mf_std_line]
                    mf_all_std_list.append(mf_std_list)
        mf_all_std_list = np.array(mf_all_std_list)
        
        higher_index = set()
        for i, row in enumerate(mf_all_std_list):
            for j, num in enumerate(row):
                if num > op_in["std_threshold"]:                    
                    higher_index.add(i)
        higher_index_list = list(higher_index)
        print("higher index list", list(cv_forces_list[higher_index_list]))
        mf_all_std_list_modified = np.delete(mf_all_std_list, higher_index_list, axis=0)
        cv_forces_list_modified = np.delete(cv_forces_list, higher_index_list, axis=0)
        assert len(mf_all_std_list_modified) == len(cv_forces_list_modified)
        
        task_path = Path("label_std")
        task_path.mkdir(exist_ok=True, parents=True)
        with set_directory(task_path):
            plt.figure(figsize=(8,6))
            for cv_index in range(cv_dim):
                plt.hist(mf_all_std_list_modified[:,cv_index], 100, label = "cv%s"%(cv_index+1))
            plt.legend()
            plt.xlabel("mean force std")
            plt.ylabel("frequency")
            plt.savefig(mf_std_fig)
             
            
        op_out = OPIO({
            "mf_std_fig":  task_path.joinpath(mf_std_fig),
            "cv_forces": list(cv_forces_list_modified)
        })
        return op_out

    