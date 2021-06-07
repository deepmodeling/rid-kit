from context import rid
import unittest, mock, os, shutil
from rid.enhcMD import run_enhc


class TestRunEnhc(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.path.exists("test_case/case_run_enhc"):
            shutil.rmtree("test_case/case_run_enhc")
        os.mkdir("test_case/case_run_enhc")
        os.chdir("test_case/case_run_enhc")
        for ii in range(2):
            os.mkdir("iter.{:06d}".format(ii))
            os.mkdir("iter.{:06d}/00.enhcMD".format(ii))
            for jj in range(2):
                os.mkdir("iter.{:06d}/00.enhcMD/{:03d}".format(ii, jj))
                if ii > 0:
                    open("iter.{:06d}/00.enhcMD/{:03d}/graph.pb".format(ii, jj), 'w').close()
        os.chdir("../..")


    def setUp(self):
        self.json_file = "benchmark_json/rid.json"
        self.machine_json = "benchmark_json/machine.json"
        self.base_dir = os.path.abspath("test_case/case_run_enhc")
        pass

    @mock.patch("rid.enhcMD.Submission")
    @mock.patch("rid.enhcMD.Task")
    def test_run_enhc_iter0(self, mock_task, mock_submission):
        run_enhc(0, self.json_file,self.machine_json, self.base_dir)
        self.assertEqual(mock_task.call_count, 4)
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
                'command':  "gmx mdrun -ntmpi 1 -nt 8 -plumed plumed.bf.dat 1> gmx_mdrun.log 2> gmx_mdrun.log",
                'task_work_path': "{:03d}".format(ii),
                'outlog': 'gmx_mdrun.log', 
                'errlog': 'gmx_mdrun.log'
            }
            for key in args_1.keys():
                self.assertEqual(set(mock_task.call_args_list[ii+2][1][key].split()), set(args_2[key].split()))
        self.assertEqual(mock_submission.call_count, 2)
        self.assertEqual(len(mock_submission.mock_calls), 4)
    
    @mock.patch("rid.enhcMD.Submission")
    @mock.patch("rid.enhcMD.Task")
    def test_run_enhc_iter1(self, mock_task, mock_submission):
        run_enhc(1, self.json_file,self.machine_json, self.base_dir)
        self.assertEqual(mock_task.call_count, 4)
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
                'command':  "gmx mdrun -ntmpi 1 -nt 8 -plumed plumed.dat 1> gmx_mdrun.log 2> gmx_mdrun.log",
                'task_work_path': "{:03d}".format(ii),
                'outlog': 'gmx_mdrun.log', 
                'errlog': 'gmx_mdrun.log'
            }
            for key in args_1.keys():
                self.assertEqual(set(mock_task.call_args_list[ii+2][1][key].split()), set(args_2[key].split()))
        self.assertEqual(mock_submission.call_count, 2)
        self.assertEqual(len(mock_submission.mock_calls), 4)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists("test_case/case_run_enhc"):
            shutil.rmtree("test_case/case_run_enhc")


if __name__ == "__main__":
    unittest.main()