from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Parameter
    )


class AdjustTrustLevel(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "trust_lvl_1": float,
                "trust_lvl_2": float,
                "init_trust_lvl_1": float,
                "init_trust_lvl_2": float,
                "numb_cluster": int,
                "numb_cluster_threshold": int,
                "adjust_amplifier": float,  
                "max_level_multiple": float  
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "adjust_trust_lvl_1": float,
                "adjust_trust_lvl_2": float
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:
        if op_in["numb_cluster"] < op_in["numb_cluster_threshold"]:
            adjust_trust_lvl_1 = op_in["trust_lvl_1"] * op_in["adjust_amplifier"]
            adjust_trust_lvl_2 = op_in["trust_lvl_2"] * op_in["adjust_amplifier"]
        else:
            adjust_trust_lvl_1 = op_in["init_trust_lvl_1"] 
            adjust_trust_lvl_2 = op_in["init_trust_lvl_2"] 
        if adjust_trust_lvl_1 > op_in["init_trust_lvl_1"] * op_in["max_level_multiple"]:
            adjust_trust_lvl_1 = op_in["init_trust_lvl_1"] 
            adjust_trust_lvl_2 = op_in["init_trust_lvl_2"] 
        op_out = OPIO({
            "adjust_trust_lvl_1": adjust_trust_lvl_1,
            "adjust_trust_lvl_2": adjust_trust_lvl_2
        })
        return op_out

    