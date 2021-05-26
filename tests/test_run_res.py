from context import rid
import unittest, mock, os, shutil
from rid.resMD import run_res

class TestRunRes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.path.exists("test_case/case_run_res"):
            shutil.rmtree("test_case/case_run_res")
        os.mkdir("test_case/case_run_res")
        os.chdir("test_case/case_run_res")
        iter_idx = 1
        os.mkdir("iter.{:06d}".format(iter_idx))
        os.mkdir("iter.{:06d}/01.resMD".format(iter_idx))
        for jj in range(2):
            os.mkdir("iter.{:06d}/01.resMD/{:03d}".format(iter_idx, jj))
        os.chdir("../..")


    def setUp(self):
        self.json_file = "benchmark_json/rid.json"
        self.base_dir = os.path.abspath("test_case/case_run_res")
        self.machine_json = "benchmark_json/machine.json"
        pass

    @mock.patch("rid.resMD.Submission")
    @mock.patch("rid.resMD.Task")
    def test_run_res_iter0(self, mock_task, mock_submission):
        run_res(1, self.json_file,self.machine_json, self.base_dir)
        # print(input("123456"))
        for ii in range(2):
            args_1 = {
                'command': "gmx grompp -maxwarn 1 1> gmx_grompp.log 2> gmx_grompp.log",
                'task_work_path': "{:03d}".format(ii),
                'outlog': 'gmx_grompp.log', 
                'errlog': 'gmx_grompp.log'
            }
            for key in args_1.keys():
                self.assertEqual(set(mock_task.call_args_list[ii][1][key].split()), set(args_1[key].split()))
        for ii in range(2):
            args_2 = {
                'command':  "gmx mdrun -ntmpi 1 -nt 8 -plumed plumed.res.dat 1> gmx_mdrun.log 2> gmx_mdrun.log",
                'task_work_path': "{:03d}".format(ii),
                'outlog': 'gmx_mdrun.log', 
                'errlog': 'gmx_mdrun.log'
            }
            for key in args_1.keys():
                self.assertEqual(set(mock_task.call_args_list[ii+2][1][key].split()), set(args_2[key].split()))
        self.assertEqual(mock_task.call_count, 4)
        self.assertEqual(mock_submission.call_count, 2)
        self.assertEqual(len(mock_submission.mock_calls), 4)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists("test_case/case_run_res"):
            shutil.rmtree("test_case/case_run_res")


if __name__ == "__main__":
    unittest.main()