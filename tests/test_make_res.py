from context import rid
import unittest, os, shutil, json, glob, filecmp
from rid.resMD import make_res
from convert import convert_pbtxt_to_pb

angle_out_content = """-2.503469 2.463779
-2.999999 3.000000
-2.397668 2.357936
"""
class TestMakeRes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.path.exists("test_case/case_make_res"):
            shutil.rmtree("test_case/case_make_res")
        cwd = os.getcwd()
        ff_list = [os.path.abspath(kk) for kk in glob.glob("benchmark_mol/*.gro")]
        graph_txt_list = [os.path.abspath(kk) for kk in glob.glob("benchmark_case/000/*.pbtxt")]
        for ii, gtf in enumerate(graph_txt_list):
            convert_pbtxt_to_pb(gtf, "benchmark_case/000/graph.{:03d}.pb".format(ii))
        graph_list = [os.path.abspath(kk) for kk in glob.glob("benchmark_case/000/*.pb")]
        for jj in range(2):
            for its in range(2):
                os.makedirs("test_case/case_make_res/test_iter{}/iter.{:06d}/00.enhcMD/{:03d}".format(jj, jj, its))
                os.chdir("test_case/case_make_res/test_iter{}/iter.{:06d}/00.enhcMD/{:03d}".format(jj, jj, its))
                os.mkdir("confs")
                for ii in range(3):
                    os.symlink(ff_list[ii], "confs/conf{}.gro".format(ii))
                with open("angle.rad.out", 'w') as aro:
                    aro.write(angle_out_content)
                os.symlink("confs/conf0.gro", "conf.gro")
                os.symlink("confs/conf0.gro", "confout.gro")
                if jj == 1:
                    for gf in range(len(graph_list)):
                        os.symlink(graph_list[gf], "graph.{:03d}.pb".format(gf))
                        if not os.path.exists("../graph.{:03d}.pb".format(gf)):
                            os.symlink(graph_list[gf], "../graph.{:03d}.pb".format(gf))
                    os.chdir("../../..")
                    if not os.path.exists("cluster_threshold.dat"):
                        with open("cluster_threshold.dat", "w") as clust:
                            clust.write("1.500000")
                os.chdir(cwd)

    def setUp(self):
        self.json = "benchmark_json/rid_make_res.json"
        self.cv_file = "benchmark_json/cv_make_res.json"
        self.mol_path = "benchmark_mol"
        self.test_dir_0 = os.path.abspath("test_case/case_make_res/test_iter0")
        self.test_dir_1 = os.path.abspath("test_case/case_make_res/test_iter1")
        self.cwd = os.getcwd()
        
    def test_make_res_iter0(self):
        make_res(0, self.json, self.cv_file, self.mol_path, "test_case/case_make_res/test_iter0")
        os.chdir(self.test_dir_0)
        self.assertTrue(os.path.exists("cluster_threshold.dat"))
        os.chdir("iter.000000/01.resMD")
        sel_dir = glob.glob("000.*")
        res_file = ["centers.out", "grompp.mdp", "plumed.res.dat", "topol.top", "conf.gro"]
        enhc_file = ["cls.sel.out", "cluster_threshold.dat", "num_of_cluster.dat", "sel.out"]
        for ii in sel_dir:
            os.chdir(ii)
            for rf in res_file:
                self.assertTrue(os.path.exists(rf))
            os.chdir("..")
        os.chdir(self.cwd)
        os.chdir(self.test_dir_0 + "/iter.000000/00.enhcMD/000")
        for ef in enhc_file:
            self.assertTrue(os.path.exists(ef))
        os.chdir(self.cwd)
    
    def test_make_res_iter1(self):
        make_res(1, self.json, self.cv_file, self.mol_path, "test_case/case_make_res/test_iter1")
        os.chdir(self.test_dir_1)
        self.assertTrue(os.path.exists("cluster_threshold.dat"))
        os.chdir("iter.000001/01.resMD")
        sel_dir = glob.glob("000.*")
        for ii in sel_dir:
            os.chdir(ii)
            self.assertTrue(os.path.exists("centers.out"))
            self.assertTrue(os.path.exists("grompp.mdp"))
            self.assertTrue(os.path.exists("plumed.res.dat"))
            self.assertTrue(os.path.exists("topol.top"))
            self.assertTrue(os.path.exists("conf.gro"))
            os.chdir("..")
        os.chdir(self.cwd)
        pass

    @classmethod
    def tearDownClass(cls):
        if os.path.exists("test_case/case_make_res"):
            shutil.rmtree("test_case/case_make_res")
        graph_list = [os.path.abspath(kk) for kk in glob.glob("benchmark_case/000/*.pb")]
        for gf in graph_list:
            os.remove(gf)
        pass
    

if __name__ == "__main__":
    unittest.main()