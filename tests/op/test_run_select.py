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

class Test_MockedRunSelect(unittest.TestCase):
    def setUp(self):
        self.taskname = "000"
        self.datapath = "data"
        cls_data = np.loadtxt(Path(self.datapath)/"cls_sel.out")
        np.save(Path(self.datapath)/"cls_sel.out.npy", cls_data)
        cls_idx = np.loadtxt(Path(self.datapath)/"cls_sel.ndx", dtype=int)
        np.save(Path(self.datapath)/"cls_sel.ndx.npy", cls_idx)
    
    def tearDown(self):
        ii = Path(self.taskname)
        if ii.is_dir():
            shutil.rmtree(ii)
        cls_out_npy = Path(self.datapath)/"cls_sel.out.npy"
        cls_ndx_npy = Path(self.datapath)/"cls_sel.ndx.npy"
        for file in [cls_ndx_npy, cls_out_npy]:
            if os.path.exists(file):
                os.remove(file)
        
    @patch('rid.op.run_select.slice_xtc')
    def test(self, mocked_run):
        mocked_run.return_value = None
        op = RunSelect()
        data = Path(self.datapath)
        cluster_indx = data/"cls_sel.ndx.npy"
        cluster_data = data/"cls_sel.out.npy"
        xtc_traj = data/"traj_comp.xtc"
        top_path = data/"topol.top"
        model_path1 = ".."/data/"models/model_000.pb"
        model_path2 = ".."/data/"models/model_001.pb"
        model_path3 = ".."/data/"models/model_002.pb"
        gmx_config = {"type":"gmx","nsteps": 50,"method":"restrained", "output_freq": 1, "temperature": 300, "kappas": [500,500],
                      "dt": 0.002, "output_mode": "both", "ntmpi": 1, "nt": 8, "max_warning": 0}

        op_in1 = OPIO(
            {
                "task_name": self.taskname,
                "cluster_selection_index": cluster_indx,
                "cluster_selection_data": cluster_data,
                "models": None,
                "trust_lvl_1": 0.02,
                "trust_lvl_2": 0.03,
                "xtc_traj": xtc_traj,
                "topology": top_path,
                "label_config": gmx_config,
                "dt": 0.002,
                "output_freq": 2500,
                "slice_mode": "gmx"
            }
        )
        op_in2 = OPIO(
            {
                "task_name": self.taskname,
                "cluster_selection_index": cluster_indx,
                "cluster_selection_data": cluster_data,
                "models": None,
                "trust_lvl_1": 0.02,
                "trust_lvl_2": 0.03,
                "xtc_traj": xtc_traj,
                "topology": top_path,
                "label_config": gmx_config,
                "dt": 0.002,
                "output_freq": 2500,
                "slice_mode": "mdtraj"
            }
        )
        op_in3 = OPIO(
            {
                "task_name": self.taskname,
                "cluster_selection_index": cluster_indx,
                "cluster_selection_data": cluster_data,
                "models": [model_path1, model_path2, model_path3],
                "trust_lvl_1": 0.02,
                "trust_lvl_2": 0.03,
                "xtc_traj": xtc_traj,
                "topology": top_path,
                "label_config": gmx_config,
                "dt": 0.002,
                "output_freq": 2500,
                "slice_mode": "gmx"
            }
        )
        op_in4 = OPIO(
            {
                "task_name": self.taskname,
                "cluster_selection_index": cluster_indx,
                "cluster_selection_data": cluster_data,
                "models": None,
                "trust_lvl_1": 0.02,
                "trust_lvl_2": 0.03,
                "xtc_traj": xtc_traj,
                "topology": top_path,
                "label_config": gmx_config,
                "dt": 0.002,
                "output_freq": 2500,
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