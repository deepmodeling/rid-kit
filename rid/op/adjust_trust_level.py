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
                "number_of_cluster": int,
                "number_of_cluster_threshold": int,
                "amplifier": float,  # 1.5 by default
                "max_level_multiple": float  # 8 by default
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
        if op_in["number_of_cluster"] < op_in["number_of_cluster_threshold"]:
            adjust_trust_lvl_1 = op_in["trust_lvl_1"] * op_in["amplifier"]
            adjust_trust_lvl_2 = op_in["trust_lvl_2"] * op_in["amplifier"]
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

    