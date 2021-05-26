from context import rid
import unittest, shutil, os
from rid.train import make_train

def create_file(file_path, content=None):
    fp = open(file_path, 'w')
    if content is not None:
        fp.write(str(content))
    fp.close()

def create_path(fpath):
    if not os.path.exists(fpath):
        os.makedirs(fpath)

class TestMakeTrain(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.path.exists("tmp_make_train"):
            shutil.rmtree("tmp_make_train")
        for ii in range(2):
            os.makedirs("tmp_make_train/iter.{:06d}/02.train".format(ii))
            os.makedirs("tmp_make_train/iter.{:06d}/01.resMD".format(ii))
            create_file("tmp_make_train/iter.{:06d}/01.resMD/data.raw".format(ii), content = ii)

    def setUp(self):
        self.json = "benchmark_json/rid.json"
        self.base_dir = "tmp_make_train"
        self.cwd = os.getcwd()
        self.ckfile = ["checkpoint", "model.ckpt.data-00000-of-00001", "model.ckpt.index", "model.ckpt.meta"]

    def test_make_train_iter0(self):
        make_train(0, self.json, self.base_dir)
        work_path = self.base_dir + "/iter.000000/02.train"
        os.chdir(work_path)
        data_file = ["data.new.raw", "data.old.raw", "data.raw"]
        self.assertTrue(all([os.path.exists("data/" + ff) for ff in data_file]))
        fc = open("data/data.raw", 'r')
        ret = fc.read()
        fc.close()
        self.assertTrue(ret == '0')
        for ii in range(4):
            self.assertTrue(os.path.exists("{:03d}".format(ii)))
            os.chdir("{:03d}".format(ii))
            self.assertTrue(all([os.path.exists("data/" + ff) for ff in data_file]))
            os.chdir("..")
        os.chdir(self.cwd)
        print("\ncurrent path\n", os.getcwd())

    def test_make_train_iter1(self):
        create_path("tmp_make_train/iter.000000/02.train/data")
        if not os.path.exists("tmp_make_train/iter.000000/02.train/data/data.raw"):
            create_file("tmp_make_train/iter.000000/02.train/data/data.raw", content="9")

        for ii in range(4):
            create_path("tmp_make_train/iter.000000/02.train/{:03d}/".format(ii))
            for ff in self.ckfile:
                create_file("tmp_make_train/iter.000000/02.train/{:03d}/{}".format(ii, ff))
        make_train(1, self.json, self.base_dir)
        work_path = self.base_dir + "/iter.000001/02.train"
        os.chdir(work_path)
        data_file = ["data.new.raw", "data.old.raw", "data.raw"]
        self.assertTrue(all([os.path.exists("data/" + ff) for ff in data_file]))
        fc = open("data/data.raw", 'r')
        ret = fc.read()
        fc.close()
        self.assertTrue(ret == '01')
        for ii in range(4):
            self.assertTrue(os.path.exists("{:03d}".format(ii)))
            os.chdir("{:03d}".format(ii))
            self.assertTrue(all([os.path.exists("data/" + ff) for ff in data_file]))
            self.assertTrue(all([os.path.exists(ff) for ff in self.ckfile]))
            os.chdir("..")
        os.chdir(self.cwd)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists("tmp_make_train"):
            shutil.rmtree("tmp_make_train")


if __name__ == "__main__":
    unittest.main()