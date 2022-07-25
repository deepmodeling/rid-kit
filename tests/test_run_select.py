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
from rid.op.run_select import RunSelect
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
        self.taskname = "selected"
        self.datapath = "data"
    
    def tearDown(self):
        ii = Path(self.taskname)
        if ii.is_dir():
            shutil.rmtree(ii)
    
    @patch('rid.op.run_select.slice_xtc')
    def test(self, mocked_run):
        mocked_run.return_value = None
        op = RunSelect()
        data = Path(self.datapath)
        cluster_indx = data/"cls_sel.ndx"
        cluster_data = data/"cls_sel.out"
        xtc_traj = data/"traj_comp.xtc"
        top_path = data/"topol.top"
        model_path1 = ".."/data/"model_000.pb"
        model_path2 = ".."/data/"model_001.pb"
        model_path3 = ".."/data/"model_002.pb"

        op_in1 = OPIO(
            {
                "task_name": self.taskname,
                "culster_selection_index": cluster_indx,
                "culster_selection_data": cluster_data,
                "models": None,
                "trust_lvl_1": 0.02,
                "trust_lvl_2": 0.03,
                "xtc_traj": xtc_traj,
                "topology": top_path,
                "dt": 0.002,
                "slice_mode": "gmx"
            }
        )
        op_in2 = OPIO(
            {
                "task_name": self.taskname,
                "culster_selection_index": cluster_indx,
                "culster_selection_data": cluster_data,
                "models": None,
                "trust_lvl_1": 0.02,
                "trust_lvl_2": 0.03,
                "xtc_traj": xtc_traj,
                "topology": top_path,
                "dt": 0.002,
                "slice_mode": "mdtraj"
            }
        )
        op_in3 = OPIO(
            {
                "task_name": self.taskname,
                "culster_selection_index": cluster_indx,
                "culster_selection_data": cluster_data,
                "models": [model_path1, model_path2, model_path3],
                "trust_lvl_1": 0.02,
                "trust_lvl_2": 0.03,
                "xtc_traj": xtc_traj,
                "topology": top_path,
                "dt": 0.002,
                "slice_mode": "gmx"
            }
        )
        op_in4 = OPIO(
            {
                "task_name": self.taskname,
                "culster_selection_index": cluster_indx,
                "culster_selection_data": cluster_data,
                "models": None,
                "trust_lvl_1": 0.02,
                "trust_lvl_2": 0.03,
                "xtc_traj": xtc_traj,
                "topology": top_path,
                "dt": 0.002,
                "slice_mode": "others"
            }
        )
    
        op_out1 = op.execute(op_in1)
        op_out2 = op.execute(op_in2)
        op_out3 = op.execute(op_in3)

        self.assertTrue(op_out1["model_devi"])
        self.assertTrue(op_out2["model_devi"])
        self.assertTrue(op_out3["model_devi"])
        self.assertTrue(op_out1["selected_indices"])
        self.assertTrue(op_out2["selected_indices"])
        self.assertTrue(op_out3["selected_indices"])
        self.assertRaises(RuntimeError, op.execute, op_in4)