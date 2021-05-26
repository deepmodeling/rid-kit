from context import rid
import unittest, mock, os, shutil, glob
from rid.train import run_train
from rid.lib.nn import NN_PATH


def create_file(file_path, content=None):
    fp = open(file_path, 'w')
    if content is not None:
        fp.write(str(content))
    fp.close()

def create_path(fpath):
    if not os.path.exists(fpath):
        os.makedirs(fpath)

class TestRunTrain(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.path.exists("tmp_run_train"):
            shutil.rmtree("tmp_run_train")
        for ii in range(2):
            os.makedirs("tmp_run_train/iter.{:06d}/02.train/data".format(ii))
            create_file("tmp_run_train/iter.{:06d}/02.train/data/data.raw".format(ii), content = ii)
            for jj in range(2):
                os.makedirs("tmp_run_train/iter.{:06d}/02.train/{:03d}".format(ii, jj))

    def setUp(self):
        self.json_file = "benchmark_json/rid_test_train.json"
        self.base_dir = os.path.abspath("tmp_run_train")
        self.machine_json = "benchmark_json/machine.json"
        self.cv_file = "benchmark_json/cv.json"
        if len(glob.glob("tmp_run_train/iter.000000/02.train/[0-9]*[0-9]/graph.pb")) < 2:
            create_file("tmp_run_train/iter.000000/02.train/000/graph.pb")
            create_file("tmp_run_train/iter.000000/02.train/001/graph.pb")

    
    @mock.patch("rid.train.cal_cv_dim", return_value=[10, 10])
    @mock.patch("rid.train.Submission")
    @mock.patch("rid.train.Task")
    def test_run_train_raw(self, mock_task, mock_submission, mock_cv):
        graph_file = glob.glob("tmp_run_train/iter.000000/02.train/*.pb")
        if len(graph_file) > 0:
            for ff in graph_file:
                os.remove(ff)
        create_file("tmp_run_train/iter.000000/02.train/data/data.new.raw", content = 0)
        run_train(0, self.json_file, self.machine_json, self.cv_file, self.base_dir)
        self.assertEqual(mock_task.call_count, 4)
        for ii in range(2):
            args_1 = {
                'command': "python3 {}/train.py -t 8 --resnet  -n 10 11 12 13  -c 10 10  -b 128 -e 2000 -l 0.0008 --decay-steps 120 --decay-rate 0.96 1> train.log 2> train.log".format(NN_PATH),
                'task_work_path': "{:03d}".format(ii),
                'outlog': 'train.log', 
                'errlog': 'train.log'
            }
            for key in args_1.keys():
                self.assertEqual(set(mock_task.call_args_list[ii][1][key].split()), set(args_1[key].split()))
        for ii in range(2):
            args_2 = {
                'command':  "python3 {}/freeze.py -o graph.pb 1> freeze.log 2> freeze.log".format(NN_PATH),
                'task_work_path': "{:03d}".format(ii),
                'outlog': 'freeze.log', 
                'errlog': 'freeze.log'
            }
            for key in args_1.keys():
                self.assertEqual(set(mock_task.call_args_list[ii+2][1][key].split()), set(args_2[key].split()))
        self.assertEqual(mock_submission.call_count, 2)
        self.assertEqual(len(mock_submission.mock_calls), 4)
        graph_file = glob.glob("tmp_run_train/iter.000000/02.train/*.pb")
        self.assertEqual(len(graph_file), 2)
    
    @mock.patch("rid.train.cal_cv_dim")
    @mock.patch("rid.train.Submission")
    @mock.patch("rid.train.Task")
    def test_run_train_new(self, mock_task, mock_submission, mock_cv):
        if len(glob.glob("tmp_run_train/iter.000000/02.train/*.pb")) < 2:
            create_file("tmp_run_train/iter.000000/02.train/graph.000.pb")
            create_file("tmp_run_train/iter.000000/02.train/graph.001.pb")
        create_file("tmp_run_train/iter.000001/02.train/data/data.new.raw", content=None)
        run_train(1, self.json_file, self.machine_json, self.cv_file, self.base_dir)
        mock_task.assert_not_called()
        mock_submission.assert_not_called()
        mock_cv.assert_not_called()
        graph_file = glob.glob("tmp_run_train/iter.000001/02.train/*.pb")
        self.assertEqual(len(graph_file), 2)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists("tmp_run_train"):
            shutil.rmtree("tmp_run_train")


if __name__ == "__main__":
    unittest.main()