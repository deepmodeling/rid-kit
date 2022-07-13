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
from rid.constants import tf_model_name
from rid.nn.train_net import train
from rid.nn.freeze import freeze_model
from rid.utils import load_txt


class TrainModel(OP):

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
        data_shape = load_txt(op_in["data"]).shape
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
