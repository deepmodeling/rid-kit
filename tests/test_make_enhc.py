from context import rid
from rid.enhcMD import make_enhc
import unittest, os, shutil, json, glob, filecmp


test_files = ['angle.rad.out', 'cls.sel.angle.0.out',
              'cls.sel.angle.out', 'cls.sel.out', 
              'cluster_threshold.dat', 'conf.gro', 
              'conf_init.gro', 'confout.gro', 'confs', 
              'gmx_grompp.err', 'grompp.mdp', 
              'num_of_cluster.dat', 'plm.out', 
              'plumed.bf.dat', 'plumed.dat', 'posre.itp', 
              'sel.angle.out', 'sel.out', 
              'topol.top', 'traj_comp.xtc', 'trust_lvl1.dat']

class TestMakeRid(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cwd = os.getcwd()
        if os.path.exists('./test_case/case_make_enhc_iter0'):
            shutil.rmtree('./test_case/case_make_enhc_iter0')
        os.mkdir('./test_case/case_make_enhc_iter0')
        
        if os.path.exists('./test_case/case_make_enhc_iter1'):
            shutil.rmtree('./test_case/case_make_enhc_iter1')
        for it in range(2):
            os.makedirs('./test_case/case_make_enhc_iter1/iter.{:06d}'.format(it))
        os.chdir('./test_case/case_make_enhc_iter1/iter.000000')
        os.makedirs("02.train")
        for ii in range(2):
            os.makedirs("00.enhcMD/{:03d}".format(ii))
            for ff in test_files:
                if ff == "confout.gro" or ff == "conf_init.gro":
                    fp = open("00.enhcMD/{:03d}/{}".format(ii, ff), "w")
                    fp.write(str(ff))
                    fp.close()
                else:
                    fp = open("00.enhcMD/{:03d}/{}".format(ii, ff), "w").close()
            with open("02.train/graph.{:03d}.pb".format(ii), 'w') as pb:
                pb.write(str(ii))
        os.chdir(cwd)


    def setUp(self):
        self.test_dir_0 = './test_case/case_make_enhc_iter0'
        self.test_dir_1 = './test_case/case_make_enhc_iter1'
        self.test_rid_out = './test_case/case_make_enhc_iter0'
        self.benchmark_json_dir = "benchmark_json"
        self.benchmark_mol_dir = "benchmark_mol"
        self.cwd = os.getcwd()
        self.conf_file_list = sorted([ os.path.join(self.cwd, pp) for pp in glob.glob(self.benchmark_mol_dir + "/conf*gro")])
        self.json_file = os.path.join(self.benchmark_json_dir, "rid.json")
        self.cv_file = os.path.join(self.benchmark_json_dir, "cv.json")
        fp = open (self.json_file, 'r')
        self.jdata = json.load (fp)
        fp.close()
        
    def test_make_enhc_iter0(self):
        idx_iter = 0
        make_enhc(
            iter_index=idx_iter, 
            json_file=self.json_file, 
            graph_files=[], 
            mol_dir=self.benchmark_mol_dir, 
            cv_file=self.cv_file, 
            base_dir=self.test_dir_0)
        enhc_files = ['topol.top','grompp.mdp','conf.gro', 'conf_init.gro', 'trust_lvl1.dat']
        enhc_files += ['plumed.bf.dat', 'plumed.dat']
        iter_dir = self.test_rid_out + "/iter.000000"
        self.assertTrue(os.path.exists(iter_dir))
        enhc_dir = iter_dir + "/00.enhcMD"
        walker_dir_list = glob.glob(enhc_dir + "/[0-9]*[0-9]")
        self.assertTrue(len(walker_dir_list) == int(self.jdata["numb_walkers"]))
        
        for idx_walker, walker_path in enumerate(walker_dir_list):
            os.chdir(walker_path)
            for enhcf in enhc_files:
                self.assertTrue(os.path.exists(enhcf))
            self.assertTrue(filecmp.cmp("conf.gro", self.conf_file_list[idx_walker]))
            self.assertTrue(filecmp.cmp("conf_init.gro", self.conf_file_list[idx_walker]))
            os.chdir(self.cwd)
        pass

    def test_make_enhc_iter1(self):
        idx_iter = 1
        prev_model = glob.glob('./test_case/case_make_enhc_iter1/iter.000000/02.train/*.pb')
        make_enhc(
            iter_index=idx_iter, 
            json_file=self.json_file, 
            graph_files=prev_model, 
            mol_dir=self.benchmark_mol_dir, 
            cv_file=self.cv_file, 
            base_dir=self.test_dir_1)
        iter_dir = self.test_rid_out + "/iter.000001"
        enhc_dir = iter_dir + "/00.enhcMD"
        walker_dir_list = glob.glob(enhc_dir + "/[0-9]*[0-9]")   
        prev_conf_out = os.path.join(self.test_dir_1, 'iter.000000/00.enhcMD/{:03d}/confout.gro')     
        for idx_walker, walker_path in enumerate(walker_dir_list):
            os.chdir(walker_path)
            self.assertTrue(filecmp.cmp("conf.gro", prev_conf_out.format(idx_walker)))
            self.assertTrue(filecmp.cmp("conf_init.gro", self.conf_file_list[idx_walker]))
            self.assertTrue(len(glob.glob("*.pb")) > 0 )
            os.chdir(self.cwd)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists('./test_case/case_make_enhc_iter0'):
            shutil.rmtree('./test_case/case_make_enhc_iter0')
        if os.path.exists('./test_case/case_make_enhc_iter1'):
            shutil.rmtree('./test_case/case_make_enhc_iter1')
        

if __name__ == "__main__":
    unittest.main()