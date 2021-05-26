from context import rid
import unittest, os, shutil, glob
import numpy as np
from rid.lib import make_ndx
from rid.lib import make_def

conf_content = """    1ASN      N    1   1.120   2.777   1.604 -0.6201  0.2837  0.1286
    1ASN     H1    2   1.019   2.783   1.601 -0.5687  1.3644  0.5561
    1ASN     H2    3   1.153   2.821   1.689  1.0401  0.0871 -0.3980
    1ASN     H3    4   1.161   2.838   1.535  1.5769  1.9049  2.8145
    1ASN     CA    5   1.168   2.638   1.579  0.2191 -0.0103  0.2400
    1ASN     HA    6   1.125   2.593   1.668 -1.8332 -1.8518 -1.6038
    1ASN     CB    7   1.106   2.584   1.452  0.2401  0.2985  0.6630
    1ASN    HB1    8   1.159   2.612   1.361 -0.9504  1.6209  0.3733
    1ASN    HB2    9   0.999   2.603   1.451 -0.2925 -2.3867 -2.0511
    1ASN     CG   10   1.113   2.433   1.472  0.0870  0.4633  0.9903
    1ASN    OD1   11   1.042   2.370   1.550  0.1301 -0.3259  0.2855
    1ASN    ND2   12   1.205   2.366   1.401  0.4309  0.6028 -0.7961
    1ASN   HD21   13   1.272   2.414   1.342  1.2710 -1.2210 -1.3453
    1ASN   HD22   14   1.222   2.269   1.424  5.5803  1.7282  0.7256
    1ASN      C   15   1.318   2.626   1.577  0.3732 -0.4514 -0.9065
    1ASN      O   16   1.362   2.526   1.628  0.5213 -0.7825  0.1721
    2LEU      N   17   1.388   2.728   1.523 -0.3547 -0.1980 -0.5603
    2LEU      H   18   1.345   2.810   1.483  0.6576  0.1328 -1.0069
    2LEU     CA   19   1.533   2.725   1.497 -0.4584 -0.3686  1.2633
    2LEU     HA   20   1.556   2.636   1.438  0.7515  2.6616 -3.0643
    2LEU     CB   21   1.566   2.829   1.389 -0.0831  0.3878  0.2438
    2LEU    HB1   22   1.487   2.817   1.315  2.2682  2.4597 -2.7386
    2LEU    HB2   23   1.659   2.802   1.340 -0.2791 -1.3193  0.7917
    2LEU     CG   24   1.568   2.974   1.437  0.3395  0.1510 -0.9432
    2LEU     HG   25   1.486   2.998   1.505  1.2657  1.2880 -0.1959
    2LEU    CD1   26   1.692   3.007   1.524  0.1025  0.2489  0.0084
    2LEU   HD11   27   1.668   2.984   1.628  0.1905 -0.2929 -0.0880
    2LEU   HD12   28   1.783   2.961   1.486 -0.3061  0.8661 -1.7743
    2LEU   HD13   29   1.707   3.114   1.516 -0.5369  0.3690  0.3310
    2LEU    CD2   30   1.569   3.074   1.323 -0.7748 -0.7214  0.3664
    2LEU   HD21   31   1.549   3.181   1.333 -1.9083 -0.7859 -1.0137
    2LEU   HD22   32   1.669   3.075   1.277  0.2388  2.3545  2.4601
    2LEU   HD23   33   1.493   3.032   1.257 -1.3688 -1.5016  1.5433
    2LEU      C   34   1.624   2.724   1.623  0.3777 -1.1878 -0.4136
    2LEU      O   35   1.742   2.689   1.611  0.6090 -0.2264  0.6847
    3SOL    HW2 9472   0.342   4.525   3.339 -1.6886  0.6562  1.4661
    4CL      CL 9473   1.227   1.871   2.756  0.0210 -0.1631 -0.3716"""


class TestGetNdx(unittest.TestCase):
    def setUp(self):
        self.mol_dir = "benchmark_mol"
        self.line = "    1ASN      N    1   1.120   2.777   1.604 -0.6201  0.2837  0.1286"
        
    def test_get_res_idx(self):
        res_idx = make_ndx.get_res_idx(self.line)
        self.assertTrue(res_idx == 1)

    def test_get_res_name(self):
        res_name = make_ndx.get_res_name(self.line)
        self.assertTrue(res_name == "ASN")

    def test_get_atom_idx(self):
        atom_idx = make_ndx.get_atom_idx(self.line)
        self.assertTrue(atom_idx == 1)

    def test_get_atom_name(self):
        atom_name = make_ndx.get_atom_name(self.line)
        self.assertTrue(atom_name == "N")

class TestMakeNdx(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.path.exists("tmp_lib_make/"):
            shutil.rmtree("tmp_lib_make")
        os.mkdir("tmp_lib_make")
        with open("tmp_lib_make/test.gro", 'w') as test_gro:
            test_gro.write("test\n4\n")
            test_gro.write(conf_content)
            test_gro.write("\n   4.57054   4.57054   4.57054")

    def setUp(self):
        self.test_dir = "tmp_lib_make"
        self.res_list = conf_content.split("\n")

    def test_make_residue_atoms(self):
        res_id = make_ndx.make_residue_atoms(self.res_list, 0, 16)
        self.assertTrue(res_id == {'N': 1, 'CA': 5, 'CB': 7, 'C': 15, 'O': 16})
        pass

    def test_make_ndx(self):
        residues, residue_atoms = make_ndx.make_ndx(self.test_dir + "/test.gro")
        benchmark_residues = [['ASN', 0, 16], ['LEU', 16, 35], ['SOL', 35, 36], ['CL', 36, 37]]
        benchmark_residue_atoms = [{'N': 1, 'CA': 5, 'CB': 7, 'C': 15, 'O': 16}, {'N': 17, 'CA': 19, 'CB': 21, 'C': 34, 'O': 35}, {}, {}]
        self.assertTrue(residues == benchmark_residues)
        self.assertTrue(residue_atoms == benchmark_residue_atoms)
    
    def test_make_protein_atom_index(self):
        tmp_index = make_ndx.make_protein_atom_index(self.test_dir + "/test.gro")
        benchmark_index = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
        self.assertTrue(tmp_index == benchmark_index)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists("tmp_lib_make/"):
            shutil.rmtree("tmp_lib_make")


class TestMakeDef(unittest.TestCase):
    def setUp(self):
        self.residues = [['ASN', 0, 16], ['LEU', 16, 35], ['SOL', 35, 36], ['CL', 36, 37]]
        self.residue_atoms = [{'N': 1, 'CA': 5, 'CB': 7, 'C': 15, 'O': 16}, {'N': 17, 'CA': 19, 'CB': 21, 'C': 34, 'O': 35}, {}, {}]

    def test_make_angle(self):
        dih_angles = [[ {'name': ['C'], 'resid_shift': -1},
                        {'name': ['N'], 'resid_shift': 0},
                        {'name': ['CA'], 'resid_shift': 0}, 
                        {'name': ['C'], 'resid_shift': 0}], 
                      [ {'name': ['N'], 'resid_shift': 0}, 
                        {'name': ['CA'], 'resid_shift': 0}, 
                        {'name': ['C'], 'resid_shift': 0}, 
                        {'name': ['N'], 'resid_shift': 1} ]]
        angle_names, angle_atom_idxes = make_def.make_general_angle_def(self.residue_atoms, dih_angles, fmt_alpha="%03d", fmt_angle="%02d")
        self.assertTrue(angle_names == ['dih-000-01', 'dih-001-00'])
        self.assertTrue(angle_atom_idxes == [[1, 5, 15, 17], [15, 17, 19, 34]])
        ret = make_def.make_angle_def (angle_names, angle_atom_idxes)
        self.assertTrue(ret == "dih-000-01: TORSION ATOMS=1,5,15,17 \ndih-001-00: TORSION ATOMS=15,17,19,34 \n")


    def test_make_dist(self):
        dist_names, dist_atom_idxes = make_def.make_general_dist_def (self.residues, self.residue_atoms,sel_residue_names=["LEU", "ASN"], sel_atom_names=["N", "CA"],fmt_residue = "%03d",exclude = 0.1)
        self.assertTrue(dist_names == ['dist-000-001'])
        self.assertTrue(dist_atom_idxes == [[1, 17]])
        ret = make_def.make_dist_def (dist_names, dist_atom_idxes)
        self.assertTrue(ret == "dist-000-001: DISTANCE ATOMS=1,17 \n")
        dist_names, dist_atom_idxes = make_def.make_general_dist_def (self.residues, self.residue_atoms,sel_residue_names=["LEU", "ASN"], sel_atom_names=["N", "CA"],fmt_residue = "%03d",exclude = 1000)
        self.assertTrue(dist_names == [])
        self.assertTrue(dist_atom_idxes == [])
        
        
if __name__ == "__main__":
    unittest.main()