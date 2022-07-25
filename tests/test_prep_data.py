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
from rid.op.prep_data import CollectData,MergeData
from rid.utils import load_txt, save_txt, set_directory
from pathlib import Path
from rid.constants import (
        data_new,
        data_old,
        data_raw,
        force_out,
        center_out_name
    )

class Test_Collectdata(unittest.TestCase):
    def setUp(self):
        self.datapath = "data"
    
    def tearDown(self):
        ii=Path(data_new)
        os.remove(ii)
    
    def test(self):
        op = CollectData()
        data = Path(self.datapath)
        op_in = OPIO(
            {
                "forces": [data/force_out],
                "centers": [data/center_out_name]
            }
        )
        op_out = op.execute(op_in)
        data1 = np.load(op_out["data_new"])
        data2 = np.load((data/data_new))
        np.testing.assert_almost_equal(data1, data2,6)
        
class Test_Mergedata(unittest.TestCase):
    def setUp(self):
        self.data = "data"
        
    def tearDown(self):
        ii=Path(data_raw)
        os.remove(ii)
    
    def test(self):
        op = MergeData()
        data = Path(self.data)
        op_in1 = OPIO(
            {
                "data_old": data/data_old,
                "data_new": data/data_new,
            }
        )
        op_out1 = op.execute(op_in1)
        data1 = np.load(op_out1["data_raw"])
        data2 = np.load((data/data_raw))
        np.testing.assert_almost_equal(data1, data2, 6)
        
        op_in2 = OPIO(
            {
                "data_old": None,
                "data_new": data/data_new,
            }
        )
        op_out2 = op.execute(op_in2)
        data1 = np.load(op_out2["data_raw"])
        data2 = np.load((data/"data_new.raw.npy"))
        np.testing.assert_almost_equal(data1, data2, 6)
        
        op_in3 = OPIO(
            {
                "data_old": data/data_old,
                "data_new": data/"data2.new.npy",
            }
        )
        op_out3 = op.execute(op_in3)
        data1 = np.load(op_out3["data_raw"])
        data2 = np.load((data/"data_old.raw.npy"))
        np.testing.assert_almost_equal(data1, data2, 6)
        