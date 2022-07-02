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
        np.savetxt(data_new, data, fmt="%.6e")
        op_out = OPIO(
            {
                "data_new": data_new
            }
        )
        return op_out


class MergeData(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "data_old": Artifact(List[Path]),
                "data_new": Artifact(List[Path])
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
        if os.stat(op_in["data_new"]).st_size == 0:
            return op_in["data_new"]
        _data_old = load_txt(op_in["data_old"])
        _data_new = load_txt(op_in["data_new"])
        data = np.concatenate((_data_old, _data_new), axis=0)
        np.savetxt(data_raw, data, fmt="%.6e")
        op_out = OPIO(
            {
                "data_raw": data
            }
        )
        return op_out
