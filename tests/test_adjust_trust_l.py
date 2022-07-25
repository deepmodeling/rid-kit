import os
import numpy as np
import unittest
from pathlib import Path
from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Parameter
    )
from context import rid
from rid.op.adjust_trust_level import AdjustTrustLevel


class Test_Adjust_trust_l(unittest.TestCase):
    def test(self):
        op = AdjustTrustLevel()
        op_in1 = OPIO(
                {
                    "trust_lvl_1": 8.0,
                    "trust_lvl_2": 10.0,
                    "init_trust_lvl_1": 2.0,
                    "init_trust_lvl_2": 3.0,
                    "numb_cluster": 6,
                    "numb_cluster_threshold": 8,
                    "adjust_amplifier": 1.5,  
                    "max_level_multiple": 3  
                })
        op_in2 = OPIO(
                {
                    "trust_lvl_1": 8.0,
                    "trust_lvl_2": 10.0,
                    "init_trust_lvl_1": 2.0,
                    "init_trust_lvl_2": 3.0,
                    "numb_cluster": 10,
                    "numb_cluster_threshold": 8,
                    "adjust_amplifier": 1.5,  
                    "max_level_multiple": 3  
                })
        op_out1 = op.execute(op_in1)
        op_out2 = op.execute(op_in2)
        self.assertEqual(op_out1["adjust_trust_lvl_1"], 2.0)
        self.assertEqual(op_out1["adjust_trust_lvl_2"], 3.0)
        self.assertEqual(op_out2["adjust_trust_lvl_1"], 2.0)
        self.assertEqual(op_out2["adjust_trust_lvl_2"], 3.0)