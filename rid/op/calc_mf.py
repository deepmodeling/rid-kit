from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact
)

import json, shutil
from typing import Tuple, List, Optional, Dict, Union
from pathlib import Path
import numpy as np
from rid.constants import (
        force_out
    )
from rid.utils import load_txt, save_txt
from rid.common.gromacs.command import get_grompp_cmd, get_mdrun_cmd


class CalcMF(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "plm_out": Artifact(Path),
                "kappa": List[float],
                "at": Artifact(Path),
                "tail": float,
                "angular_mask": Optional[Union[np.ndarray, List]],
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "forces": Artifact(Path)
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:
        data = load_txt(op_in["plm_out"])
        data = data[:, 1:]  # removr the first column(time index).
        centers = load_txt(op_in["at"])

        nframes = data.shape[0]
        
        angular_boolean = (np.array(op_in["angular_mask"], dtype=int) == 1)
        init_angle = data[0][angular_boolean]
        for ii in range(1, nframes):
            current_angle = data[ii][angular_boolean]
            angular_diff = current_angle - init_angle
            current_angle[angular_diff < -np.pi] += 2 * np.pi
            current_angle[angular_diff >= np.pi] -= 2 * np.pi
            data[ii][angular_boolean] = current_angle

        start_f = int(nframes * (1-op_in["tail"]))
        avgins = np.average(data[start_f:, :], axis=0)

        # diff = np.zeros(avgins.shape)
        diff = avgins - centers
        angular_diff = diff[angular_boolean]
        angular_diff[angular_diff < -np.pi] += 2 * np.pi
        angular_diff[angular_diff >  np.pi] -= 2 * np.pi
        diff[angular_boolean] = angular_diff
        ff = np.multiply(op_in["kappa"], diff)
        np.savetxt(force_out,  np.reshape(ff, [1, -1]), fmt='%.10e')
        op_out = OPIO(
            {
                "forces": Path(force_out)
            }
        )
        return op_out
