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
from utils import DATA_NEW, DATA_OLD, DATA_EMPTY, DATA_RAW
from rid.op.prep_data import CollectData,MergeData
from pathlib import Path
from rid.constants import (
        data_new,
        data_old,
        data_raw,
        cv_force_out
    )


class Test_Collectdata(unittest.TestCase):
    def setUp(self):
        self.datapath = "data"
        np.save(Path(self.datapath)/data_new, DATA_NEW)

    def tearDown(self):
        ii=Path(data_new)
        os.remove(ii)
        os.remove(Path(self.datapath)/data_new)
    
    def test(self):
        op = CollectData()
        data = Path(self.datapath)
        op_in = OPIO(
            {
                "cv_forces": [data/cv_force_out]
            }
        )
        op_out = op.execute(op_in)
        data1 = np.load(op_out["data_new"])
        data2 = np.load((data/data_new))
        np.testing.assert_almost_equal(data1, data2, 6)
        

class Test_Mergedata(unittest.TestCase):
    def setUp(self):
        self.data = "data"
        np.save(Path(self.data)/data_new, DATA_NEW)
        np.save(Path(self.data)/data_old, DATA_OLD)
        np.save(Path(self.data)/"data2.new.npy", DATA_EMPTY)
        
    def tearDown(self):
        ii=Path(data_raw)
        os.remove(ii)
        os.remove(Path(self.data)/data_new)
        os.remove(Path(self.data)/data_old)
        os.remove(Path(self.data)/"data2.new.npy")
    
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
        np.testing.assert_almost_equal(data1, DATA_RAW, 6)
        op_in2 = OPIO(
            {
                "data_old": None,
                "data_new": data/data_new,
            }
        )
        op_out2 = op.execute(op_in2)
        data1 = np.load(op_out2["data_raw"])
        np.testing.assert_almost_equal(data1, DATA_NEW, 6)
        
        op_in3 = OPIO(
            {
                "data_old": data/data_old,
                "data_new": data/"data2.new.npy",
            }
        )
        op_out3 = op.execute(op_in3)
        data1 = np.load(op_out3["data_raw"])
        np.testing.assert_almost_equal(data1, DATA_OLD, 6)
        