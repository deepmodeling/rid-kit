from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact
)
import os
from typing import List
from pathlib import Path
from rid.constants import (
        data_new,
        data_raw
    )
from rid.utils import load_txt
import numpy as np


class CollectData(OP):

    r"""Gather data of different simulations to a single file.
    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "cv_forces": Artifact(List[Path])
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
        cv_forces = []
        for idx in range(len(op_in["cv_forces"])):
            if op_in["cv_forces"][idx]:
                cv_force = load_txt(op_in["cv_forces"][idx])
                cv_forces = np.append(cv_forces, cv_force)
        if op_in["cv_forces"]:
            cv_forces = np.reshape(cv_forces, [-1, len(cv_force)])
            data = cv_forces
        else:
            data = np.array([])
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
        _data_new = np.load(op_in["data_new"])
        if len(_data_new) == 0:
            return OPIO({"data_raw": op_in["data_old"]})
        _data_old = np.load(op_in["data_old"])
        data = np.concatenate((_data_old, _data_new), axis=0)
        np.save(data_raw, data)
        op_out = OPIO(
            {
                "data_raw": Path(data_raw)
            }
        )
        return op_out
