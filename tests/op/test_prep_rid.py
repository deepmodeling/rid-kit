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
from rid.op.prep_rid import PrepRiD
from rid.utils import load_txt, save_txt, set_directory
from pathlib import Path
import shutil
from rid.constants import (
        init_conf_gmx_name
    )

class Test_PrepRid(unittest.TestCase):
    def setUp(self):
        self.datapath = "data"
        self.confs = "conf.gro"
        self.numb_walkers = 2
    
    def tearDown(self):
        data = Path(self.datapath)
        ii=data/self.confs
        for idx in range(self.numb_walkers):
            os.remove(init_conf_gmx_name.format(idx=idx))
    
    def test(self):
        op = PrepRiD()
        data = Path(self.datapath)
        rid_path = data/"rid.json"
        conf_path = data/"conf.gro"
        op_in1 = OPIO(
            {
                "confs": [conf_path],
                "rid_config": rid_path
            }
        )
        op_in2 = OPIO(
            {
                "confs": [conf_path, conf_path],
                "rid_config": rid_path
            }
        )
        op_in3 = OPIO(
            {
                "confs": [conf_path, conf_path, conf_path],
                "rid_config": rid_path
            }
        )
        op_out1 = op.execute(op_in1)
        op_out2 = op.execute(op_in2)
        op_out3 = op.execute(op_in3)
        self.assertEqual(op_out1["numb_walkers"],2)
        self.assertEqual(len(op_out1["confs"]),2)
        self.assertEqual(len(op_out2["confs"]),2)
        self.assertEqual(len(op_out3["confs"]),2)