from context import rid
import unittest, os, shutil, glob
import numpy as np
from rid.lib.gen import gen_plumed
from rid.lib.gen import gen_mdp

benchmark_general_plumed = """WHOLEMOLECULES ENTITY0=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22

dih-001-00: TORSION ATOMS=5,7,9,15 
dih-001-01: TORSION ATOMS=7,9,15,17 
dist-001-002: DISTANCE ATOMS=7,17 

dpfe: DEEPFE TRUST_LVL_1=1.0 TRUST_LVL_2=2.0 MODEL=graph.pb ARG=dih-001-00,dih-001-01,dist-001-002 

PRINT STRIDE=5 ARG=dih-001-00,dih-001-01,dist-001-002 FILE=plm.out 


"""

class TestMakeNdx(unittest.TestCase):
    def setUp(self):
        self.conf = "benchmark_mol/conf0.gro"
        self.cv_file = "benchmark_json/cv.json"
 
    def test_general_plumed(self):
        ret = gen_plumed.general_plumed (TASK='dpbias',CONF=self.conf,JSON=self.cv_file,kappa = 500.0,temp = 3000.0, tau = 10.0, gamma = 0.1, pstride = 5, pfile = "plm.out")
        ret = "".join([kk for kk in ret.strip().split() if kk != ""])
        benchmark = [kk for kk in benchmark_general_plumed.strip().split() if kk != ""]
        for part in benchmark:
            self.assertTrue(part in ret)
        benchmark = sorted(list("".join(benchmark)))
        ret = sorted(list(ret))
        self.assertTrue(ret, benchmark)
        
        
if __name__ == "__main__":
    unittest.main()