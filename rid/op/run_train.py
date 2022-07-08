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


class TrainModel(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "task_id": int,
                "cv_dim": int,
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
        train_config = op_in["train_config"]
        train(
            cv_dim=op_in["cv_dim"],
            neurons=train_config["neurons"],
            angular_mask=op_in["angular_mask"],
            numb_threads=train_config["numb_threads"],
            resnet=train_config["resnet"],
            use_mix=train_config["use_mix"],
            restart=train_config["restart"],
            batch_size=train_config["batch_size"],
            epoches=train_config["epoches"],
            lr=train_config["learning_rate"],
            decay_steps=train_config["decay_steps"],
            decay_rate=train_config["decay_rate"],
            old_ratio=train_config["old_ratio"],
            decay_steps_inner=train_config["decay_steps_inner"],
            drop_out_rate=train_config["drop_out_rate"]
        )
        out_put_name = tf_model_name.format(idx=op_in["task_id"])
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


class __TrainModel(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "task_id": int,
                "cv_dim": int,
                "angular_mask": List,
                "data": Artifact(Path),
                "neurons": List,
                "numb_threads": Parameter(int, default=8),
                "resnet": Parameter(bool, default=True),
                "use_mix": Parameter(bool, default=False),
                "restart": Parameter(bool, default=False),
                "batch_size": Parameter(int, default=128),
                "epoches": Parameter(int, default=12000),
                "learning_rate": Parameter(float, default=0.0008),
                "decay_steps": Parameter(int, default=120),
                "decay_rate": Parameter(float, default=0.96),
                "old_ratio": Parameter(float, default=7.0),
                "decay_steps_inner": Parameter(int, default=0),
                "drop_out_rate": Parameter(float, default=0.3)
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
        
        train(
            cv_dim=op_in["cv_dim"],
            neurons=op_in["neurons"],
            angular_mask=op_in["angular_mask"],
            numb_threads=op_in["numb_threads"],
            resnet=op_in["resnet"],
            use_mix=op_in["use_mix"],
            restart=op_in["restart"],
            batch_size=op_in["batch_size"],
            epoches=op_in["epoches"],
            lr=op_in["learning_rate"],
            decay_steps=op_in["decay_steps"],
            decay_rate=op_in["decay_rate"],
            old_ratio=op_in["old_ratio"],
            decay_steps_inner=op_in["decay_steps_inner"],
            drop_out_rate=op_in["drop_out_rate"]
        )
        out_put_name = tf_model_name.format(idx=op_in["task_id"])
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