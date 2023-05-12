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
from rid.op.prep_select import PrepSelect
from rid.utils import load_txt, save_txt, set_directory
from pathlib import Path
import shutil

class Test_PrepSelect(unittest.TestCase):
    def setUp(self):
        self.taskname = "000"
        self.datapath = "data"
    
    def tearDown(self):
        ii = Path(self.taskname)
        if ii.is_dir():
            shutil.rmtree(ii)
    
    def test(self):
        op = PrepSelect()
        data = Path(self.datapath)
        plm_path = data/"plm.out"

        op_in1 = OPIO(
            {
                "task_name": self.taskname,
                "plm_out": plm_path,
                "cluster_threshold": 0.05,
                "angular_mask": [1,1],
                "weights": [1,1],
                "numb_cluster_upper": 5,
                "numb_cluster_lower": 3,
                "max_selection": 5,
                "if_make_threshold": True
            }
        )
        op_in2 = OPIO(
            {
                "task_name": self.taskname,
                "plm_out": plm_path,
                "cluster_threshold": 0.05,
                "angular_mask": [1,1],
                "weights": [1,1],
                "numb_cluster_upper": 5,
                "numb_cluster_lower": 3,
                "max_selection": 5,
                "if_make_threshold": False
            }
        )
    
        op_out1 = op.execute(op_in1)
        op_out2 = op.execute(op_in2)

        self.assertEqual(op_out1["cluster_threshold"], 0.05)
        self.assertEqual(op_out1["numb_cluster"], 5)
        self.assertEqual(op_out2["cluster_threshold"], 0.05)