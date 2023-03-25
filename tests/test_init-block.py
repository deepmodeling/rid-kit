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
    data_raw,
    walker_tag_fmt,
    tf_model_name,
    sel_ndx_name,
    gmx_conf_out,
    gmx_xtc_name,
    gmx_mdrun_log,
    model_tag_fmt
)

from rid.utils import normalize_resources

from rid.superop.exploration import Exploration
from rid.superop.label import Label
from rid.superop.selector import Selector
from rid.superop.data import DataGenerator
from rid.superop.blocks import InitBlock

from context import (
    skip_ut_with_dflow,
    skip_ut_with_dflow_reason,
    default_image,
    default_host,
)

from mocked_ops import (
    clear_files,
    clear_dir,
    make_mocked_init_models,
    make_mocked_init_confs,
    MockedPrepExplore,
    MockedRunExplore,
    MockedCheckLabelInputs,
    MockedPrepLabel,
    MockedRunLabel,
    MockedPrepSelect,
    MockedRunSelect,
    MockedCollectData,
    MockedMergeData,
    MockedTrain
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

exploration_op = Exploration(
    "exploration",
    MockedPrepExplore,
    MockedRunExplore,
    default_config,
    default_config)

label_op = Label(
    "label",
    MockedCheckLabelInputs,
    MockedPrepLabel,
    MockedRunLabel,
    default_config,
    default_config)

select_op = Selector(
    "select",
    MockedPrepSelect,
    MockedRunSelect,
    default_config,
    default_config)

data_op = DataGenerator(
    "gen-data",
    MockedCollectData,
    MockedMergeData,
    default_config)

init_block_op = InitBlock(
    "init-block",
    exploration_op,
    select_op,
    label_op,
    data_op,
    MockedTrain,
    default_config,
)
    
@unittest.skipIf(skip_ut_with_dflow, skip_ut_with_dflow_reason)
class TestMockedInitBlock(unittest.TestCase):
    def setUp(self):
        self.init_models = make_mocked_init_models(4)
        self.init_confs = make_mocked_init_confs(3)
        trust_lvl_1 = 2.0
        trust_lvl_2 = 4.0
        self.exploration_config = {"type":"gmx"}
        self.label_config = {"type":"gmx"}
        self.cv_config = {"cv_dim": 3}
        self.train_config = {"numb_models": 2}
        walker_tags = []
        model_tags = []
        self.numb_walkers = 2
        for ii in range(self.numb_walkers):
            walker_tags.append(walker_tag_fmt.format(idx=ii))
        for ii in range(self.train_config["numb_models"]):
            model_tags.append(model_tag_fmt.format(idx=ii))
        self.models = upload_artifact(self.init_models)
        self.confs = upload_artifact(self.init_confs)
        self.trust_lvl_1 = [trust_lvl_1]*self.numb_walkers
        self.trust_lvl_2 = [trust_lvl_2]*self.numb_walkers
        self.cluster_threshold = [2.0]*self.numb_walkers
        self.angular_mask = [0.]*3
        self.weights = [1.]*3
        self.numb_cluster_upper = 8
        self.numb_cluster_lower = 4
        self.max_selection = 10
        self.dt = 0.002
        self.tail = 0.9
        self.output_freq = 50
        self.slice_mode = "gmx"
        self.if_make_threshold = False
        self.data_out = "init-block"
        self.block_tag = "iter-001"
        self.walker_tags = walker_tags
        self.model_tags = model_tags

    def tearDown(self):
        clear_files(self.init_models)
        clear_files(self.init_confs)
        clear_dir(self.data_out)
            
    def test_exploration(self):
        init_block_step = Step(
            'init-block-step',
            template = init_block_op,
            parameters = {
                "block_tag" : self.block_tag,
                "walker_tags" : self.walker_tags,
                "model_tags": self.model_tags,
                "trust_lvl_1" : self.trust_lvl_1,
                "trust_lvl_2": self.trust_lvl_2,
                "exploration_config" : self.exploration_config,
                "cv_config" : self.cv_config,
                "cluster_threshold": self.cluster_threshold,
                "angular_mask": self.angular_mask,
                "weights": self.weights,
                "numb_cluster_upper": self.numb_cluster_upper,
                "numb_cluster_lower": self.numb_cluster_lower,
                "max_selection": self.max_selection,
                "dt": self.dt,
                "output_freq": self.output_freq,
                "slice_mode": self.slice_mode,
                "label_config": self.label_config,
                "tail": self.tail,
                "train_config": self.train_config
            },
            artifacts = {
                "models" : None,
                "forcefield": None,
                "topology" : None,
                "inputfile": None,
                "confs" : self.confs,
                "index_file": None,
                "dp_files": None,
                "cv_file": None
            },
        )
        
        wf = Workflow(name="init-block-test", host=default_host)
        wf.add(init_block_step)
        wf.submit()
        
        while wf.query_status() in ["Pending", "Running"]:
            time.sleep(4)

        self.assertEqual(wf.query_status(), "Succeeded")
        step = wf.query_step(name="init-block-step")[0]
        self.assertEqual(step.phase, "Succeeded")

        if not os.path.exists(self.data_out):
            os.mkdir(self.data_out)
        download_artifact(step.outputs.artifacts["conf_outs"],path=self.data_out)
        download_artifact(step.outputs.artifacts["data"],path=self.data_out)
        download_artifact(step.outputs.artifacts["exploration_md_log"],path=self.data_out)
        download_artifact(step.outputs.artifacts["exploration_trajectory"],path=self.data_out)
        download_artifact(step.outputs.artifacts["selection_index"],path=self.data_out)
        download_artifact(step.outputs.artifacts["models"],path=self.data_out)
        sub_path = "/000/"
        self.assertTrue(os.path.isfile(self.data_out+sub_path+gmx_conf_out))
        self.assertTrue(os.path.isfile(self.data_out+sub_path+gmx_mdrun_log))
        self.assertTrue(os.path.isfile(self.data_out+sub_path+gmx_xtc_name))
        self.assertTrue(os.path.isfile(self.data_out+sub_path+sel_ndx_name))
        self.assertTrue(os.path.isfile(self.data_out+"/"+data_raw))
        self.assertTrue(os.path.isfile(self.data_out+sub_path+tf_model_name.format(tag="000")))
        with open(self.data_out+sub_path+gmx_conf_out,"r") as f:
            l1 = f.readline()
            self.assertEqual(l1, "This is init conf 0")
