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
from rid.op.prep_exploration import PrepExplore
from rid.utils import load_txt, save_txt, set_directory
from pathlib import Path
import shutil
from rid.constants import (
        explore_task_pattern, 
        gmx_conf_name,
        gmx_top_name,
        gmx_mdp_name, 
        plumed_input_name,
        plumed_output_name
    )

class Test_PrepExplore(unittest.TestCase):
    def setUp(self):
        self.datapath = "data"
        self.taskname = "explored"
    
    def tearDown(self):
        ii=Path(self.taskname)
        shutil.rmtree(ii)
    
    def test(self):
        op = PrepExplore()
        data = Path(self.datapath)
        gmx_config = {"type":"gmx","nsteps": 50, "output_freq": 1, "temperature": 300, 
                      "dt": 0.002, "output_mode": "both", "ntmpi": 1, "nt": 8, "max_warning": 0}
        cv_config = {"mode": "torsion", "selected_resid": [1, 2], "cv_file": []}
        cv_config2 = {"mode": "custom", "selected_resid": [1, 2], "cv_file": ["data/plumed.dat"]}
        op_in1 = OPIO(
            {
                "models": [data/"model_000.pb", data/"model_001.pb", data/"model_002.pb"],
                "topology": data/"topol.top",
                "conf": data/"conf.gro",
                "cv_file": None,
                "trust_lvl_1": 2.0,
                "trust_lvl_2": 3.0,
                "exploration_config": gmx_config,
                "cv_config": cv_config,
                "task_name": "explored",
                "block_tag": "iter-001"
            }
        )
        op_in2 = OPIO(
            {
                "models": None,
                "topology": data/"topol.top",
                "conf": data/"conf.gro",
                "cv_file": None,
                "trust_lvl_1": 2.0,
                "trust_lvl_2": 3.0,
                "exploration_config": gmx_config,
                "cv_config": cv_config,
                "task_name": "explored",
                "block_tag":"iter-001"
            }
        )
        op_in3 = OPIO(
            {
                "models": [data/"model_000.pb", data/"model_001.pb", data/"model_002.pb"],
                "topology": data/"topol.top",
                "conf": data/"conf.gro",
                "cv_file": [data/"plumed.dat"],
                "trust_lvl_1": 2.0,
                "trust_lvl_2": 3.0,
                "exploration_config": gmx_config,
                "cv_config": cv_config2,
                "task_name": "explored",
                "block_tag":"iter-001"
            }
        )
        op_out1 = op.execute(op_in1)
        op_out2 = op.execute(op_in2)
        op_out3 = op.execute(op_in3)
        self.assertEqual(op_out1["cv_dim"],2)
        self.assertTrue(op_out1["task_path"].is_dir())
        assert gmx_conf_name in os.listdir(op_out1["task_path"]), "configuration file not in outdir"
        assert gmx_top_name in os.listdir(op_out1["task_path"]), "topology file not in outdir"
        assert plumed_input_name in os.listdir(op_out1["task_path"]), "plumed input file not in outdir"
        assert gmx_mdp_name in os.listdir(op_out1["task_path"]), "mdp file not in outdir"
        
        self.assertEqual(op_out2["cv_dim"],2)
        self.assertTrue(op_out2["task_path"].is_dir())
        assert gmx_conf_name in os.listdir(op_out2["task_path"]), "configuration file not in outdir"
        assert gmx_top_name in os.listdir(op_out2["task_path"]), "topology file not in outdir"
        assert plumed_input_name in os.listdir(op_out2["task_path"]), "plumed input file not in outdir"
        assert gmx_mdp_name in os.listdir(op_out2["task_path"]), "mdp file not in outdir"
        
        self.assertEqual(op_out3["cv_dim"],2)
        self.assertTrue(op_out3["task_path"].is_dir())
        assert gmx_conf_name in os.listdir(op_out3["task_path"]), "configuration file not in outdir"
        assert gmx_top_name in os.listdir(op_out3["task_path"]), "topology file not in outdir"
        assert plumed_input_name in os.listdir(op_out3["task_path"]), "plumed input file not in outdir"
        assert gmx_mdp_name in os.listdir(op_out3["task_path"]), "mdp file not in outdir"