from context import rid
import unittest, mock, os, shutil
from rid.resMD import post_res
import numpy as np
from rid.lib import LIB_PATH

def create_file(file_path, content=None):
    fp = open(file_path, 'w')
    if content is not None:
        fp.write(str(content))
    fp.close()

class TestPostRes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cwd = os.getcwd()
        if os.path.exists("tmp_post_res"):
            shutil.rmtree("tmp_post_res")
        os.makedirs("tmp_post_res/iter.000000/01.resMD")
        os.makedirs("tmp_post_res/iter.000001/01.resMD")

    def setUp(self):
        self.json_file = "benchmark_json/rid.json"
        self.base_dir = os.path.abspath("tmp_post_res")
        self.cv_file = "benchmark_json/cv_post_res.json"
        self.machine = "benchmark_json/machine.json"
        self.cwd = os.getcwd()
        pass
    
    @mock.patch("rid.resMD.cal_cv_dim", return_value=[2,0])
    @mock.patch("rid.resMD.Submission")
    @mock.patch("rid.resMD.Task")
    def test_post_res_iter0(self, mock_task, mock_submission, mock_cv):
        for ii in range(2):
            os.mkdir("tmp_post_res/iter.000000/01.resMD/{:03d}".format(ii))
            create_file("tmp_post_res/iter.000000/01.resMD/{:03d}/centers.out".format(ii), content="2.946328 -1.677528")
            create_file("tmp_post_res/iter.000000/01.resMD/{:03d}/force.out".format(ii), content="-7.6350055432e+00 -5.7704534368e+00")

        post_res(0, self.json_file, self.machine, self.cv_file, base_dir=self.base_dir)
        for ii in range(2):
            args_1 = {
                'command': "python3 {}/cmpf.py -c 2".format(LIB_PATH),
                'task_work_path': "{:03d}".format(ii),
                'outlog': 'cmpf.log', 
                'errlog': 'cmpf.log'
            }
            for key in args_1.keys():
                self.assertEqual(set(mock_task.call_args_list[ii][1][key].split()), set(args_1[key].split()))
        self.assertTrue(os.path.exists("tmp_post_res/iter.000000/01.resMD/data.raw"))
        self.assertEqual(mock_task.call_count, 2)
        self.assertEqual(mock_submission.call_count, 1)
        self.assertEqual(len(mock_submission.mock_calls), 2)

    @mock.patch("rid.resMD.cal_cv_dim")
    @mock.patch("rid.resMD.Submission")
    @mock.patch("rid.resMD.Task")
    def test_post_res_iter1(self, mock_task, mock_submission, mock_cv):
        post_res(1, self.json_file, self.machine, self.cv_file, base_dir=self.base_dir)
        self.assertTrue(os.path.exists("tmp_post_res/iter.000001/01.resMD/data.raw"))
        mock_task.assert_not_called()
        mock_submission.assert_not_called()
        mock_cv.assert_not_called()

    @classmethod
    def tearDownClass(cls):
        if os.path.exists("tmp_post_res"):
            shutil.rmtree("tmp_post_res")


if __name__ == "__main__":
    unittest.main()