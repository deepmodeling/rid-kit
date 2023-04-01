import numpy as np
from typing import List, Dict
from pathlib import Path
from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Parameter,
    BigParameter
)
from rid.constants import tf_model_name, train_fig, train_log
from rid.nn.train_net import train
from rid.nn.freeze import freeze_model
from matplotlib import pyplot as plt
from rid.utils import set_directory


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
                "train_config": BigParameter(Dict)
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "model": Artifact(Path),
                "train_log": Artifact(Path),
                "train_fig": Artifact(Path)
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
        task_path = Path(op_in["model_tag"])
        task_path.mkdir(exist_ok=True, parents=True)
        train_log_name = train_log.format(tag=op_in["model_tag"])
        with set_directory(task_path):
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
                data_path=str(op_in["data"]),
                log_name = train_log_name
            )
            out_put_name = tf_model_name.format(tag=op_in["model_tag"])
            train_fig_name = train_fig.format(tag=op_in["model_tag"])
            # plot loglog loss png
            loss_list = []
            epoch_list = []
            with open(train_log_name, "r") as f:
                while True:
                    line = f.readline()
                    if "running time" in line:
                        break
                    if "rid.nn.model" in line:
                        data = line.split(" ")
                        if len(loss_list) == 0:
                            epoch_list.append(int(data[10][:-1]))
                            loss_list.append(float(data[14][:-1]))
                        else:
                            epoch_list.append(int(data[8][:-1]))
                            loss_list.append(float(data[12][:-1]))

            plt.figure(figsize=(10, 8), dpi=100)
            plt.loglog(epoch_list,loss_list)
            plt.xlabel("log of training epoches")
            plt.ylabel("log of relative error")
            plt.title("loglog fig of training")
            plt.savefig(train_fig_name)
            
            freeze_model(
                model_folder=".",
                output=out_put_name
            )
        op_out = OPIO(
            {
                "model": task_path.joinpath(out_put_name),
                "train_log": task_path.joinpath(train_log_name),
                "train_fig": task_path.joinpath(train_fig_name)
            }
        )
        return op_out
