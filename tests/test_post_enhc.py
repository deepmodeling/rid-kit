from context import rid
import unittest, mock, os, shutil
from rid.enhcMD import post_enhc
import numpy as np

pesudo_content = """#! SET min_dih-128-00 -pi
#! SET max_dih-128-00 pi
 0.000000 -3.003329
 6.000000 3.041876
"""

def create_file(file_path, content=None):
    fp = open(file_path, 'w')
    if content is not None:
        fp.write(str(content))
    fp.close()

class TestRunEnhc(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cwd = os.getcwd()
        if os.path.exists("tmp_post_enhc"):
            shutil.rmtree("tmp_post_enhc")
        os.makedirs("tmp_post_enhc/iter.000000/00.enhcMD")
        os.chdir("tmp_post_enhc/iter.000000/00.enhcMD")
        for jj in range(2):
            os.mkdir("{:03d}".format(jj))
            create_file("{:03d}/plm.out".format(jj), content=pesudo_content)
        os.chdir(cwd)

    def setUp(self):
        self.json_file = "benchmark_json/rid.json"
        self.base_dir = os.path.abspath("tmp_post_enhc")
        self.machine = "benchmark_json/machine.json"
        pass

    @mock.patch("rid.enhcMD.Submission")
    @mock.patch("rid.enhcMD.Task")
    def test_post_enhc_iter0(self, mock_task, mock_submission):
        post_enhc(0, self.json_file, self.machine, self.base_dir)
        for ii in range(2):
            args_1 = {
                'command': "echo 0 | gmx trjconv -sep -f traj.trr -o confs/conf.gro -vel 1> gmx_split.log 2> gmx_split.log",
                'task_work_path': "{:03d}".format(ii),
                'outlog': 'gmx_split.log', 
                'errlog': 'gmx_split.log'
            }
            for key in args_1.keys():
                self.assertEqual(set(mock_task.call_args_list[ii][1][key].split()), set(args_1[key].split()))
            self.assertTrue(os.path.exists("tmp_post_enhc/iter.000000/00.enhcMD/{:03d}/angle.rad.out".format(ii)))
            result = np.loadtxt("tmp_post_enhc/iter.000000/00.enhcMD/{:03d}/angle.rad.out".format(ii))
            benchmark = np.array([-3.003329, 3.041876])
            for idx in range(len(benchmark)):
                self.assertEqual(result[idx], benchmark[idx])
        self.assertEqual(mock_task.call_count, 2)
        self.assertEqual(mock_submission.call_count, 1)
        self.assertEqual(len(mock_submission.mock_calls), 2)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists("tmp_post_enhc"):
            shutil.rmtree("tmp_post_enhc")


if __name__ == "__main__":
    unittest.main()