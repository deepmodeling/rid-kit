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
    cv_force_out,
    center_out_name,
    data_old
)

from rid.utils import normalize_resources

from rid.superop.data import DataGenerator

from context import (
    skip_ut_with_dflow,
    skip_ut_with_dflow_reason,
    default_image,
    default_host,
)

from mocked_ops import (
    clear_dir,
    make_mocked_init_data,
    MockedCollectData,
    MockedMergeData
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

# class TestMockedCollectData(unittest.TestCase):
#     def setUp(self):
#         self.forces = make_mocked_init_data(force_out,1)
#         self.centers = make_mocked_init_data(center_out_name,1)
#         self.expected_out = data_new

#     def tearDown(self):
#         init_data = self.init_data
#         out = self.expected_out
#         if Path(out).exists():
#             os.remove(out)
#         if Path(init_data).exists():
#             shutil.rmtree(init_data)

#     def test(self):
#         prep = MockedCollectData()
#         ip = OPIO(            {
#                 "forces": [self.forces],
#                 "centers": [self.centers]
#             })
#         op = prep.execute(ip)
#         self.assertEqual(Path(self.expected_out), op["data_new"])

# class TestMockedMergeData(unittest.TestCase):
#     def setUp(self):
#         self.init_data = make_mocked_init_data()
#         self.data_old = (self.init_data)/"a"
#         self.data_new = (self.init_data)/"b"
#         self.expected_out = data_raw

#     def tearDown(self):
#         init_data = self.init_data
#         out = self.expected_out
#         if Path(out).exists():
#             os.remove(out)
#         if Path(init_data).exists():
#             shutil.rmtree(init_data)

#     def test(self):
#         merge = MockedMergeData()
#         ip = OPIO(            {
#                 "data_old": self.data_old,
#                 "data_new": self.data_new
#             })
#         op = merge.execute(ip)
#         self.assertEqual(Path(self.expected_out), op["data_raw"])
        
@unittest.skipIf(skip_ut_with_dflow, skip_ut_with_dflow_reason)
class TestMockedData(unittest.TestCase):
    def setUp(self):
        self.cv_forces = make_mocked_init_data(cv_force_out,1)
        self.data_old = make_mocked_init_data(data_old,1)
        self.block_tag = "000"
        self.data_out = "data_generate"
        self.expected_out = data_raw

    def tearDown(self):
        clear_dir(self.cv_forces[0])
        clear_dir(self.data_old[0])
        clear_dir(self.data_out)
            
    def test_data(self):
        steps = DataGenerator(
            "data-steps",
            MockedCollectData,
            MockedMergeData,
            run_config = default_config
        )
        data_step = Step(
            'data-step',
            template = steps,
            parameters = {
                "block_tag" : self.block_tag
            },
            artifacts = {
            "cv_forces": upload_artifact(self.cv_forces[1]),
            "data_old": upload_artifact(self.data_old[1])
            },
        )
        
        wf = Workflow(name="data-generate", host=default_host)
        wf.add(data_step)
        wf.submit()
        
        while wf.query_status() in ["Pending", "Running"]:
            time.sleep(4)

        self.assertEqual(wf.query_status(), "Succeeded")
        step = wf.query_step(name="data-step")[0]
        self.assertEqual(step.phase, "Succeeded")

        if not os.path.exists(self.data_out):
            os.mkdir(self.data_out)
        download_artifact(step.outputs.artifacts["data"],path = self.data_out)
        self.assertTrue(os.path.isfile(self.data_out+"/"+self.expected_out))
        with open(self.data_out+"/"+self.expected_out,"r") as f:
            l1 = f.readline()
            self.assertEqual(l1, "this is %s 0\n"%data_old)
            l2 = f.readline()
            self.assertEqual(l2, "this is %s 0"%cv_force_out)