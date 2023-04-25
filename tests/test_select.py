import os
import unittest
import time

from dflow import (
    Workflow,
    Step,
    upload_artifact,
    download_artifact
)


try:
    from context import rid
except ModuleNotFoundError:
    # case of upload everything to argo, no context needed
    pass

from rid.constants import (
    cluster_selection_index_name,
    model_devi_name,
    sel_gro_name,
    cv_init_label,
    sel_ndx_name
)

from rid.utils import normalize_resources

from rid.superop.selector import Selector

from context import (
    skip_ut_with_dflow,
    skip_ut_with_dflow_reason,
    default_image,
    default_host,
)

from mocked_ops import (
    clear_dir,
    make_mocked_init_data,
    MockedPrepSelect,
    MockedRunSelect
)

from dflow.python import upload_packages
from rid import SRC_ROOT
upload_packages.append(SRC_ROOT)

default_config = normalize_resources({
        "template_config" : {
            "image": default_image,
            "image_pull_policy" : "IfNotPresent"
        }
    })
        
@unittest.skipIf(skip_ut_with_dflow, skip_ut_with_dflow_reason)
class TestMockedSelect(unittest.TestCase):
    def setUp(self):
        self.plm_out = make_mocked_init_data(data_name="plm_out",data_num=3)
        self.xtc_traj = make_mocked_init_data(data_name="xtc_traj",data_num=3)
        self.topology = make_mocked_init_data(data_name="topology",data_num=3)
        self.trust_lvl_1 = [2.0]*3
        self.trust_lvl_2 = [3.0]*3
        self.cluster_threshold = [2.0]*3
        self.angular_mask = [0,0,0]
        self.weights = [1,1,1]
        self.numb_cluster_upper = 8
        self.numb_cluster_lower = 4
        self.max_selection = 10
        self.dt = 0.002
        self.output_freq = 50
        self.slice_mode = "gmx"
        self.if_make_threshold = False
        self.task_name = ["000","001","002"]
        self.block_tag = "000"
        self.data_out = "select"
        self.label_config = {"type":"gmx"}

    def tearDown(self):
        clear_dir(self.plm_out[0])
        clear_dir(self.xtc_traj[0])
        clear_dir(self.topology[0])
        clear_dir(self.data_out)
            
    def test_select(self):
        steps = Selector(
            "select-steps",
            MockedPrepSelect,
            MockedRunSelect,
            prep_config= default_config,
            run_config = default_config
        )
        select_step = Step(
            'select-step',
            template = steps,
            parameters = {
                "label_config": self.label_config,
                "trust_lvl_1" : self.trust_lvl_1,
                "trust_lvl_2": self.trust_lvl_2,
                "cluster_threshold": self.cluster_threshold,
                "angular_mask": self.angular_mask,
                "weights": self.weights,
                "numb_cluster_upper": self.numb_cluster_upper,
                "numb_cluster_lower": self.numb_cluster_lower,
                "max_selection": self.max_selection,
                "dt": self.dt,
                "output_freq": self.output_freq,
                "slice_mode": self.slice_mode,
                "if_make_threshold": self.if_make_threshold,
                "task_names" : self.task_name,
                "block_tag" : self.block_tag
            },
            artifacts = {
                "models" : None,
                "plm_out": upload_artifact(self.plm_out[1]),
                "xtc_traj": upload_artifact(self.xtc_traj[1]),
                "topology": upload_artifact(self.topology[1])
            },
        )
        
        wf = Workflow(name="select-test", host=default_host)
        wf.add(select_step)
        wf.submit()
        
        while wf.query_status() in ["Pending", "Running"]:
            time.sleep(4)

        self.assertEqual(wf.query_status(), "Succeeded")
        step = wf.query_step(name="select-step")[0]
        self.assertEqual(step.phase, "Succeeded")
        
        if not os.path.exists(self.data_out):
            os.mkdir(self.data_out)
        download_artifact(step.outputs.artifacts["selected_confs"],path=self.data_out)
        download_artifact(step.outputs.artifacts["model_devi"],path=self.data_out)
        download_artifact(step.outputs.artifacts["selected_cv_init"],path=self.data_out)
        download_artifact(step.outputs.artifacts["selected_indices"],path=self.data_out)
        download_artifact(step.outputs.artifacts["cluster_selection_index"],path=self.data_out)
        sub_path = "/000/"
        self.assertTrue(os.path.isfile(self.data_out+sub_path+sel_gro_name.format(walker=0,idx=0)))
        self.assertTrue(os.path.isfile(self.data_out+sub_path+"cls_"+model_devi_name))
        self.assertTrue(os.path.isfile(self.data_out+sub_path+cv_init_label.format(walker=0,idx=0)))
        self.assertTrue(os.path.isfile(self.data_out+sub_path+sel_ndx_name))
        self.assertTrue(os.path.isfile(self.data_out+sub_path+cluster_selection_index_name))
        with open(self.data_out+sub_path+sel_gro_name.format(walker=0,idx=0),"r") as f:
            l1 = f.readline()
            self.assertEqual(l1, "This is mocked conf_000_0.gro")