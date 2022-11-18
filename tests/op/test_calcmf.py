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
from rid.op.calc_mf import CalcMF
from rid.utils import load_txt, save_txt, set_directory
from pathlib import Path
import shutil
from rid.constants import (
        force_out,
        plumed_output_name
    )

class Test_Calcmf(unittest.TestCase):
    def setUp(self):
        self.task_name = "force"
        self.datapath = "data"
    
    def tearDown(self):
        ii = self.task_name
        ii=Path(ii)
        if ii.is_dir():
            shutil.rmtree(ii)
    
    def test(self):
        op = CalcMF()
        data = Path(self.datapath)
        op_in = OPIO(
            {
                "task_name": "force",
                "plm_out": data/plumed_output_name,
                "kappas": [500,500],
                "at": data/"centers.out",
                "tail": 0.9,
                "angular_mask": [1, 1],
            }
        )
        op_out = op.execute(op_in)
        force1 = [float(i) for i in op_out["forces"].read_text().split(" ")]
        force2 = [float(i) for i in (data/force_out).read_text().split(" ")]
        self.assertListEqual(force1, force2)