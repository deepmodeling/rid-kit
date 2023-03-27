import os
import numpy as np
import unittest
from pathlib import Path
from dflow.python import (
    OPIO
    )
from context import rid
from rid.op.run_train import TrainModel
from rid.constants import data_raw
from pathlib import Path
from utils import DATA_RAW
import shutil

class Test_RunTrain(unittest.TestCase):
    def setUp(self):
        self.datapath = "data"
        np.save(Path(self.datapath)/data_raw, DATA_RAW)

    
    def tearDown(self):
        os.remove(Path(self.datapath)/data_raw)
        ii = Path("000")
        if ii.is_dir():
            shutil.rmtree(ii)
    
    def test(self):
        op = TrainModel()
        data = Path(os.path.abspath(self.datapath))
        train_config = {"neurons": [20, 20], "resnet": True, "batch_size": 2,
                        "epoches": 200, "init_lr": 0.0008, "decay_steps": 120,
                        "decay_rate": 0.96, "train_thread": 8, "drop_out_rate": 0.3, 
                        "numb_threads": 8, "use_mix": False, "restart": False, "decay_steps_inner": 120}
        op_in1 = OPIO(
            {
                "model_tag": "000",
                "angular_mask": [1,1],
                "data": data/data_raw,
                "train_config": train_config
            }
        )
        op_out1 = op.execute(op_in1)
        self.assertTrue(op_out1["model"])