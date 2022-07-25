from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact
)
import os
import json, shutil
from typing import Tuple, List, Optional, Dict
from pathlib import Path
from rid.constants import (
        data_new,
        data_raw
    )
from rid.task.builder import EnhcMDTaskBuilder
from rid.utils import load_txt, save_txt
import numpy as np


class CollectData(OP):

    r"""Gather data of different simulations to a single file.
    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "forces": Artifact(List[Path]),
                "centers": Artifact(List[Path])
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "data_new": Artifact(Path),
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
        ) -> OPIO:
        forces = []
        centers = []
        for idx in range(len(op_in["forces"])):
            force = load_txt(op_in["forces"][idx])
            center = load_txt(op_in["centers"][idx])
            forces = np.append(forces, force)
            centers = np.append(centers, center)
        forces = np.reshape(forces, [len(op_in["forces"]), -1])
        centers = np.reshape(centers, [len(op_in["forces"]), -1])
        data = np.concatenate((centers, forces), axis=1)
        np.save(data_new, data)
        op_out = OPIO(
            {
                "data_new": Path(data_new)
            }
        )
        return op_out


class MergeData(OP):
    
    r"""Merge old data and new generated data. 
    If old data not existed, it will return new data.
    If new data is empty, it will return old data.
    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "data_old": Artifact(Path, optional=True),
                "data_new": Artifact(Path),
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "data_raw": Artifact(Path),
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
        ) -> OPIO:
        if op_in["data_old"] is None:
            return OPIO({"data_raw": op_in["data_new"]})
        if os.stat(op_in["data_new"]).st_size == 0:
            return OPIO({"data_raw": op_in["data_old"]})
        _data_old = np.load(op_in["data_old"])
        _data_new = np.load(op_in["data_new"])
        data = np.concatenate((_data_old, _data_new), axis=0)
        np.save(data_raw, data)
        op_out = OPIO(
            {
                "data_raw": Path(data_raw)
            }
        )
        return op_out
