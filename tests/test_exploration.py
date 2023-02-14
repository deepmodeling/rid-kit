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
    walker_tag_fmt,
    plumed_output_name,
    gmx_conf_out,
    gmx_xtc_name,
    gmx_mdrun_log
)

from rid.utils import normalize_resources

from rid.superop.exploration import Exploration

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
    MockedRunExplore
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

# class TestMockedPrepExplore(unittest.TestCase):
#     def setUp(self):
#         self.models = make_mocked_init_models(4)
#         self.init_data = make_mocked_init_data()
#         self.confs = make_mocked_init_confs(1)
#         self.topology = self.init_data/"a"
#         self.cv_file = [self.init_data/"b"]
#         self.trust_lvl_1 = 2
#         self.trust_lvl_2 = 4
#         self.exploration_config = {"type":"gmx"}
#         self.cv_config = {"dim":2}
#         self.task_name = "explored"
        
#     def tearDown(self):
#         clear_data(self.models,self.init_data,self.confs,self.task_name)

#     def test(self):
#         prep = MockedPrepExplore()
#         ip = OPIO({
#                 "models": self.models,
#                 "topology": self.topology,
#                 "conf": self.confs[0],
#                 "cv_file": self.cv_file,
#                 "trust_lvl_1": self.trust_lvl_1,
#                 "trust_lvl_2": self.trust_lvl_2,
#                 "exploration_config": self.exploration_config,
#                 "cv_config": self.cv_config,
#                 "task_name": self.task_name
#                 })
#         op = prep.execute(ip)
#         self.assertTrue(Path(self.task_name))

# class TestMockedRunExplore(unittest.TestCase):
#     def setUp(self):
#         self.taskname = "explored"
#         self.task_path = Path(self.taskname)
#         self.exploration_config = {"type":"gmx"}
#         if not os.path.exists(self.taskname):
#             os.makedirs(self.taskname)
#         with set_directory(self.task_path):
#             with open(gmx_conf_name, "w") as f:
#                 f.write("conf 0")
    
#     def tearDown(self):
#         if os.path.exists(self.task_path):
#             shutil.rmtree(self.task_path)
    
#     def test(self):
#         op = MockedRunExplore()
#         op_in = OPIO(
#             {
#                 "task_path": self.task_path,
#                 "exploration_config": self.exploration_config,
#                 "models": None,
#                 "forcefield": None,
#                 "index_file": None,
#                 "dp_files": None,
#                 "cv_file": None,
#                 "inputfile": None
#             }
#             )
#         op_out = op.execute(op_in)
#         self.assertTrue(op_out["plm_out"])
#         self.assertTrue(op_out["md_log"])
#         self.assertTrue(op_out["trajectory"])
#         self.assertTrue(op_out["conf_out"])

@unittest.skipIf(skip_ut_with_dflow, skip_ut_with_dflow_reason)
class TestMockedExploration(unittest.TestCase):
    def setUp(self):
        self.init_models = make_mocked_init_models(4)
        trust_lvl_1 = 2.0
        trust_lvl_2 = 4.0
        exploration_config = {"type":"gmx"}
        cv_config = {"cv_dim": 3}
        walker_tags = []
        self.numb_walkers = 3
        self.init_confs = make_mocked_init_confs(self.numb_walkers)
        for ii in range(self.numb_walkers):
            walker_tags.append(walker_tag_fmt.format(idx=ii))
        self.models = upload_artifact(self.init_models)
        self.confs = upload_artifact(self.init_confs)
        self.trust_lvl_1 = [trust_lvl_1]*self.numb_walkers
        self.trust_lvl_2 = [trust_lvl_2]*self.numb_walkers
        self.exploration_config = exploration_config
        self.cv_config = cv_config
        self.data_out = "exploration"
        self.block_tag = "000"
        self.task_name = walker_tags

    def tearDown(self):
        clear_files(self.init_models)
        clear_files(self.init_confs)
        clear_dir(self.data_out)
            
    def test_exploration(self):
        steps = Exploration(
            "exploration",
            MockedPrepExplore,
            MockedRunExplore,
            prep_config = default_config,
            run_config = default_config
        )
        exploration_step = Step(
            'exploration-step',
            template = steps,
            parameters = {
                "trust_lvl_1" : self.trust_lvl_1,
                "trust_lvl_2": self.trust_lvl_2,
                "exploration_config" : self.exploration_config,
                "cv_config" : self.cv_config,
                "task_names" : self.task_name,
                "block_tag" : self.block_tag
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
        
        wf = Workflow(name="exploration-test", host=default_host)
        wf.add(exploration_step)
        wf.submit()
        
        while wf.query_status() in ["Pending", "Running"]:
            time.sleep(4)

        self.assertEqual(wf.query_status(), "Succeeded")
        step = wf.query_step(name="exploration-step")[0]
        self.assertEqual(step.phase, "Succeeded")

        if not os.path.exists(self.data_out):
            os.mkdir(self.data_out)
        download_artifact(step.outputs.artifacts["conf_outs"],path=self.data_out)
        download_artifact(step.outputs.artifacts["plm_out"],path=self.data_out)
        download_artifact(step.outputs.artifacts["md_log"],path=self.data_out)
        download_artifact(step.outputs.artifacts["trajectory"],path=self.data_out)
        sub_path = "/000/"
        self.assertTrue(os.path.isfile(self.data_out+sub_path+gmx_conf_out))
        self.assertTrue(os.path.isfile(self.data_out+sub_path+plumed_output_name))
        self.assertTrue(os.path.isfile(self.data_out+sub_path+gmx_mdrun_log))
        self.assertTrue(os.path.isfile(self.data_out+sub_path+gmx_xtc_name))
        with open(self.data_out+sub_path+gmx_conf_out,"r") as f:
            l1 = f.readline()
            self.assertEqual(l1, "This is init conf 0")