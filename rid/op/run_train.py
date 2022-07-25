from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Parameter
)

import numpy as np
import json, shutil
from typing import Tuple, List, Optional, Dict
from pathlib import Path
from rid.constants import tf_model_name
from rid.nn.train_net import train
from rid.nn.freeze import freeze_model
from rid.utils import load_txt


class TrainModel(OP):

    """`TrainModel` trains a set of neural network models (set by `numb_model` in `train_config`). 
    RiD-kit is powered by TensorFlow framework. The output model files are frozen in `.pb` formats by `rid.nn.freeze`.
    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "model_tag": str,
                "angular_mask": List,
                "data": Artifact(Path),
                "train_config": Dict
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "model": Artifact(Path)
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

            - `model_tag`: (`str`) Tags for neural network model files. In formats of `model_{model_tag}.pb`.
            - `angular_mask`: (`List`) Angular mask for periodic collective variables. 1 represents periodic, 0 represents non-periodic.
            - `data`: (`Artifact(Path)`) Data files for training. Prepared by `rid.op.prep_data`.
                `data` has the shape of `[number_conf, 2 * dimension_cv]` and contains the CV values and corresponding mean forces.
            - `train_config`: (`Dict`) Configuration to train neural networks, including training strategy and network structures.
          
        Returns
        -------
            Output dict with components:
        
            - `model`: (`Artifact(Path)`) Neural network models in `.pb` formats.
        """

        data_shape = np.load(op_in["data"]).shape
        cv_dim = int(data_shape[1] // 2)
        train_config = op_in["train_config"]
        train(
            cv_dim=cv_dim,
            neurons=train_config["neurons"],
            angular_mask=op_in["angular_mask"],
            numb_threads=train_config.get("numb_threads", 8),
            resnet=train_config["resnet"],
            use_mix=train_config["use_mix"],
            restart=train_config.get("restart", False),
            batch_size=train_config["batch_size"],
            epoches=train_config["epoches"],
            lr=train_config["init_lr"],
            decay_steps=train_config["decay_steps"],
            decay_rate=train_config["decay_rate"],
            drop_out_rate=train_config["drop_out_rate"],
            data_path=str(op_in["data"])
        )
        out_put_name = tf_model_name.format(tag=op_in["model_tag"])
        freeze_model(
            model_folder=".",
            output=out_put_name
        )
        op_out = OPIO(
            {
                "model": Path(out_put_name)
            }
        )
        return op_out
