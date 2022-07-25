import os
import numpy as np
import unittest
from pathlib import Path
from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Parameter
    )
from context import rid
from rid.op.run_train import TrainModel
from rid.utils import load_txt, save_txt, set_directory
from pathlib import Path
import shutil
from rid.constants import (
        explore_task_pattern, 
        gmx_conf_name,
        gmx_top_name,
        gmx_mdp_name, 
        plumed_input_name,
        plumed_output_name,
        sel_gro_name,
        init_conf_name
    )

class Test_RunTrain(unittest.TestCase):
    def setUp(self):
        self.datapath = "data"
    
    def tearDown(self):
        os.remove("./model_000.pb")
        os.remove("./model.ckpt.index")
        os.remove("./model.ckpt.meta")
        os.remove("./model.ckpt.data-00000-of-00001")
        os.remove("./checkpoint")
    
    def test(self):
        op = TrainModel()
        data = Path(self.datapath)
        train_config = {"neurons": [20, 20], "resnet": True, "batch_size": 2,
                        "epoches": 200, "init_lr": 0.0008, "decay_steps": 120,
                        "decay_rate": 0.96, "train_thread": 8, "drop_out_rate": 0.3, 
                        "numb_threads": 8, "use_mix": False, "restart": False, "decay_steps_inner": 120}
        op_in1 = OPIO(
            {
                "model_tag": "000",
                "angular_mask": [1,1],
                "data": data/"data.raw.npy",
                "train_config": train_config
            }
        )
    
        op_out1 = op.execute(op_in1)

        self.assertTrue(op_out1["model"])