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
from rid.op.prep_label import CheckLabelInputs, PrepLabel
from rid.utils import load_txt, save_txt, set_directory
from pathlib import Path
import shutil
from rid.constants import (
        plumed_input_name
    )
class Test_CheckLabelInputs(unittest.TestCase):
    def setUp(self) -> None:
        self.datapath = "data"
    
    def test(self):
        op = CheckLabelInputs()
        data = Path(self.datapath)
        conf_path = data/"conf.gro"
        conf_tags_path = data/"conf.json"
        
        op_in1 = OPIO(
            {
                "confs": [conf_path],
                "conf_tags": [conf_tags_path]
            }
        )
        op_in2 = OPIO(
            {
                "confs": None
            }
        )
        op_out1 = op.execute(op_in1)
        op_out2 = op.execute(op_in2)
        self.assertTrue(op_out1["if_continue"] == 1)
        self.assertTrue(op_out2["if_continue"] == 0)


class Test_PrepLabel(unittest.TestCase):
    def setUp(self):
        self.taskname = "labeled"
        self.datapath = "data"
    
    def tearDown(self):
        ii = Path(self.taskname)
        if ii.is_dir():
            shutil.rmtree(ii)
    
    def test(self):
        op = PrepLabel()
        data = Path(self.datapath)
        conf_path = data/"conf.gro"
        top_path = data/"topol.top"
        center_path = data/"centers.out"
        gmx_config1 = {"type":"gmx","nsteps": 50,"method":"restrained", "output_freq": 1, "temperature": 300, "kappas": [500,500],
                    "dt": 0.002, "output_mode": "both", "ntmpi": 1, "nt": 8, "max_warning": 0}
        gmx_config2 = {"type":"gmx","nsteps": 50,"method":"constrained", "output_freq": 1, "temperature": 300, "kappas": [500,500],
                    "dt": 0.002, "output_mode": "both", "ntmpi": 1, "nt": 8, "max_warning": 0}
        cv_config1 = {"mode": "torsion", "selected_resid": [1, 2], "angular_mask": [1,1],"cv_file": [""]}
        cv_config2 = {"mode": "distance", "selected_atomid": [[1,2],[3,4]], "angular_mask": [0,0],"cv_file": [""]}

        op_in1 = OPIO(
            {
                "topology": top_path,
                "conf": conf_path,
                "label_config": gmx_config1,
                "cv_file": None,
                "cv_config": cv_config1,
                "task_name": self.taskname,
                "at": center_path
            }
        )
        op_in2 = OPIO(
            {
                "topology": top_path,
                "conf": conf_path,
                "label_config": gmx_config2,
                "cv_file": None,
                "cv_config": cv_config2,
                "task_name": self.taskname
            }
        )
    
        op_out1 = op.execute(op_in1)

        with open(self.taskname+"/"+plumed_input_name, "r") as f:
            while True:
                l1 = f.readline()
                if l1 is not None:
                    if l1[0] == 'r':
                        break
            l2 = f.readline()
        a1 = float(l1.split("=")[-1].strip())
        a2 = float(l2.split("=")[-1].strip())
        c1 = np.loadtxt(center_path)[0]
        c2 = np.loadtxt(center_path)[1]
        self.assertEqual(a1, c1)
        self.assertEqual(a2, c2)
        
        op_out2 = op.execute(op_in2)
        
        with open(self.taskname+"/"+plumed_input_name, "r") as f:
            while True:
                l1 = f.readline()
                if l1 is not None:
                    if l1[0] == 'd':
                        break
            l2 = f.readline()
        a1 = l1.split("=")[-1].strip()
        a2 = l2.split("=")[-1].strip()
        c1 = "1,2"
        c2 = "3,4"
        self.assertEqual(a1, c1)
        self.assertEqual(a2, c2)