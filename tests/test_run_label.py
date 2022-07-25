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

class Test_MockedRunSelect(unittest.TestCase):
    def setUp(self):
        self.taskname = "001"
    
    def tearDown(self):
        ii = Path(self.taskname)
        if ii.is_dir():
            shutil.rmtree(ii)
    
    @patch('rid.op.run_label.run_command')
    def test(self, mocked_run):
        mocked_run.return_value = 0, 0, 1
        op = RunLabel()
        op_in1 = OPIO(
            {
                "task_path": Path(self.taskname),
                "gmx_config": {"nsteps": 50, "output_freq": 1, "temperature": 300, 
                               "dt": 0.002, "output_mode": "both", "ntmpi": 1, "nt": 8, "max_warning": 0}
            }
        )
    
        op_out1 = op.execute(op_in1)

        self.assertTrue(op_out1["plm_out"])
        self.assertTrue(op_out1["md_log"])