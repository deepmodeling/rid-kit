import os
import numpy as np
import unittest
from mock import mock, patch
from pathlib import Path
from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Parameter
    )
from context import rid
from rid.op.run_label import RunLabel
from rid.utils import load_txt, save_txt, set_directory
from pathlib import Path
import shutil

class Test_MockedRunLabel(unittest.TestCase):
    def setUp(self):
        self.taskname = "001"
        self.datapath = "data"
        os.mkdir(self.taskname)
        plm_data = np.loadtxt(Path(self.datapath)/"plm_label.out")
        np.savetxt(Path(self.taskname)/"plm.out", plm_data)
    
    def tearDown(self):
        ii = Path(self.taskname)
        if ii.is_dir():
            shutil.rmtree(ii)
    
    @patch('rid.op.run_label.run_command')
    def test(self, mocked_run):
        mocked_run.return_value = 0, 0, 1
        op = RunLabel()
        gmx_config = {"type":"gmx","nsteps": 50,"method":"restrained", "output_freq": 1, "temperature": 300, "kappas": [500,500],
                      "dt": 0.002, "output_mode": "both", "ntmpi": 1, "nt": 8, "max_warning": 0}
        cv_config = {"mode": "torsion", "selected_resid": [1, 2],"angular_mask": [ 1, 1 ],"weights": [ 1, 1 ],"cv_file": []}
        op_in1 = OPIO(
            {
                "task_path": Path(self.taskname),
                "label_config": gmx_config,
                "cv_config": cv_config,
                "task_name": self.taskname,
                "at": None,
                "forcefield": None,
                "index_file": None,
                "dp_files": None,
                "cv_file": None,
                "inputfile": None
            }
        )
    
        op_out1 = op.execute(op_in1)

        self.assertTrue(op_out1["plm_out"])
        self.assertTrue(op_out1["md_log"])