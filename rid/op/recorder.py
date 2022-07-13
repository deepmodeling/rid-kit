from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Parameter
    )
from typing import Optional
from rid.constants import block_tag_fmt

class Recorder(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "iteration": Parameter(Optional[int], default=None),
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "next_iteration": int,
                "block_tag": str
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:
        if op_in["iteration"] is None:
            next_iteration = 0
        else:
            next_iteration = op_in["iteration"] + 1
        block_tag = block_tag_fmt.format(idx_iter=next_iteration)
        op_out = OPIO({
            "next_iteration": next_iteration,
            "block_tag": block_tag
        })
        return op_out

    