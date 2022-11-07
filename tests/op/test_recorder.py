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
from rid.op.recorder import Recorder
from rid.constants import block_tag_fmt


class Test_Recorder(unittest.TestCase):
    def test(self):
        op = Recorder()
        op_in1 = OPIO(
                {
                    "iteration": None
                })
        op_in2 = OPIO(
                {
                    "iteration": 1
                })
        op_out1 = op.execute(op_in1)
        op_out2 = op.execute(op_in2)
        self.assertEqual(op_out1["next_iteration"], 1)
        self.assertEqual(op_out1["block_tag"], block_tag_fmt.format(idx_iter=1))
        self.assertEqual(op_out2["next_iteration"], 2)
        self.assertEqual(op_out2["block_tag"], block_tag_fmt.format(idx_iter=2))