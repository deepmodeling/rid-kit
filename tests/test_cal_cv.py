from context import rid
import unittest
from rid.lib.cal_cv_dim import cal_cv_dim


class TestCalCv(unittest.TestCase):
    def setUp(self):
        self.json = "benchmark_json/rid.json"
        self.cv_file = "benchmark_json/cv_make_res.json"
        self.cv_file_2 = "benchmark_json/cv.json"
        self.mol_path = "benchmark_mol"
        pass

        
    def test_cal_cv(self):
        dim_list = cal_cv_dim(self.mol_path + "/conf0.gro", self.cv_file)
        self.assertTrue(dim_list[0] == 2)
        self.assertTrue(dim_list[1] == 0)
        dim_list = cal_cv_dim(self.mol_path + "/conf0.gro", self.cv_file_2)
        self.assertTrue(dim_list[0] == 2)
        self.assertTrue(dim_list[1] == 1)


if __name__ == "__main__":
    unittest.main()