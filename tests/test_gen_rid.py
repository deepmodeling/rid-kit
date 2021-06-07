from context import rid
import unittest, os, shutil


class TestMakeRid(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.path.exists("tmp_rid"):
            shutil.rmtree("tmp_rid")
        os.mkdir("tmp_rid/")

    def setUp(self):
        self.test_dir = "tmp_rid"
        self.test_rid_out = "tmp_rid/test_out"
        self.benchmark_json_dir = "benchmark_json"
        self.benchmark_mol_dir = "benchmark_mol"
        self.json_file = os.path.join(self.benchmark_json_dir, "rid.json")
        # self.benchmark_dir = "benchmark_rid"
    
    def test_gen_rid(self):
        rid.gen_rid(out_dir=self.test_rid_out, mol_dir=self.benchmark_mol_dir, rid_json=self.json_file)
        self.assertTrue(os.path.exists(self.test_rid_out))

    def test_gen_rid_existed(self):
        self.test_rid_out_exists = "tmp_rid/test_existed"
        if not os.path.exists(self.test_rid_out_exists):
            os.mkdir(self.test_rid_out_exists)
        fp = open(self.test_rid_out_exists + "/existed.txt", 'w').close()
        rid.gen_rid(out_dir=self.test_rid_out_exists, mol_dir=self.benchmark_mol_dir, rid_json=self.json_file)
        self.assertTrue(os.path.exists(self.test_rid_out))
        self.assertTrue(os.path.exists(self.test_rid_out_exists + ".bk000/existed.txt"))

    def test_gen_rid_not_enough_conf(self):
        self.test_rid_not_enough = "tmp_rid/test_ne"
        json_not_enough = os.path.join(self.benchmark_json_dir, "rid_not_enough.json")
        with self.assertRaises(AssertionError):
            rid.gen_rid(out_dir=self.test_rid_not_enough, mol_dir=self.benchmark_mol_dir, rid_json=json_not_enough)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree("tmp_rid")

if __name__ == "__main__":
    unittest.main()