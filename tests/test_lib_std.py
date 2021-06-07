from context import rid
import unittest, os, shutil, glob
import numpy as np
from rid.lib import std
from rid.lib import LIB_PATH
from convert import convert_pbtxt_to_pb


class TestStd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cwd = os.getcwd()
        if os.path.exists("test_case/case_std"):
            shutil.rmtree("test_case/case_std")
        os.makedirs("test_case/case_std/000")
        graph_txt_list = [os.path.abspath(kk) for kk in glob.glob("benchmark_case/000/*.pbtxt")]
        for ii, gtf in enumerate(graph_txt_list):
            convert_pbtxt_to_pb(gtf, "benchmark_case/000/graph.{:03d}.pb".format(ii))
        ff_list = [os.path.abspath(kk) for kk in glob.glob("benchmark_case/000/*")]
        os.chdir("test_case/case_std/000")
        for ff in ff_list:
            os.symlink(ff, os.path.basename(ff))
        pbtxt_list = glob.glob("*.pbtxt")
        for pbtxt in pbtxt_list:
            os.remove(pbtxt)
        os.chdir(cwd)

    def setUp(self):
        self.case_dir = "test_case/case_std"
        self.walker_dir = os.path.join(self.case_dir, "000")
        self.cv_dim = 2
        self.cwd = os.getcwd()
        if os.path.exists(os.path.join(self.walker_dir, "sel.out")):
            os.remove(os.path.join(self.walker_dir, "sel.out"))
        if os.path.exists(os.path.join(self.walker_dir, "sel.angle.out")):
            os.remove(os.path.join(self.walker_dir, "sel.angle.out"))
        
    def test_std(self):
        os.chdir(self.walker_dir)
        model_list = glob.glob("*.pb")
        std.make_std(cv_dim=self.cv_dim, dataset="angle.rad.out", models=model_list, threshold=2.0, output="sel.out", output_angle="sel.angle.out")
        angle_sel = np.loadtxt("sel.angle.out")
        idx_sel = np.loadtxt("sel.out")
        benchmark_angle = np.loadtxt("benchmark.sel.angle.out")
        benchmark_idx = np.loadtxt("benchmark.sel.out")
        os.chdir(self.cwd)
        self.assertTrue(np.all(angle_sel == benchmark_angle))
        self.assertTrue(np.all(idx_sel == benchmark_idx))

    def test_std_cmd(self):
        os.chdir(self.walker_dir)
        os.system("python3 {}/std.py -m *.pb -d angle.rad.out -t 2.0 -o sel.out --output-angle sel.angle.out -c {}".format(LIB_PATH, self.cv_dim))
        angle_sel = np.loadtxt("sel.angle.out")
        idx_sel = np.loadtxt("sel.out")
        benchmark_angle = np.loadtxt("benchmark.sel.angle.out")
        benchmark_idx = np.loadtxt("benchmark.sel.out")
        os.chdir(self.cwd)
        self.assertTrue(np.all(angle_sel == benchmark_angle))
        self.assertTrue(np.all(idx_sel == benchmark_idx))
        pass
    
    def tearDown(self):
        if os.path.exists(os.path.join(self.walker_dir, "sel.out")):
            os.remove(os.path.join(self.walker_dir, "sel.out"))
        if os.path.exists(os.path.join(self.walker_dir, "sel.angle.out")):
            os.remove(os.path.join(self.walker_dir, "sel.angle.out"))
    
    @classmethod
    def tearDownClass(cls):
        if os.path.exists("test_case/case_std"):
            shutil.rmtree("test_case/case_std")
        graph_list = [os.path.abspath(kk) for kk in glob.glob("benchmark_case/000/*.pb")]
        for gf in graph_list:
            os.remove(gf)


if __name__ == "__main__":
    unittest.main()