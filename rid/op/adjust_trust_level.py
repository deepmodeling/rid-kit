from dflow.python import (
    OP,
    OPIO,
    OPIOSign
    )


class AdjustTrustLevel(OP):

    r"""AdjustTrustLeve OP adjust trust level according to the number of cluster (`numb_cluster`) and 
    the number of cluster threshold (`numb_cluster_threshold`). If numb_cluster < numb_cluster_threshold,
    the trust level will be increased. 
    Trust levels won't increase infinitly. When current trust_lvl_1 > `max_level_multiple` * init_trust_lvl_1, 
    where init_trust_lvl_1 is the initial trust level, trust_lvl_1 will be tuned to its initial value.

    """

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
        r"""Execute the OP.

        Parameters
        ----------
        op_in : dict
            Input dict with components:
        
            - `trust_lvl_1`: (`float`) Trust level 1 in the current iteration, or e0.
            - `trust_lvl_2`: (`float`) Trust level 2 in the current iteration, or e1.
            - `init_trust_lvl_1`: (`float`) Initial value of trust level 1.
            - `init_trust_lvl_1`: (`float`) Initial value of trust level 2.
            - `numb_cluster`: (`int`) Number of clusters got from cluster op.
            - `numb_cluster_threshold`: (`int`) Threshold for cluster number to adjust trust level adaptively.
                if `numb_cluster` > `numb_cluster_threshold`, `trust_lvl` will be increased.
            - `adjust_amplifier`: (`float`) Increasing multiple for trust level.
            - `max_level_multiple`: (`float`) The max multiple than trust level can be increased.

        Returns
        -------
            Output dict with components:
        
            - `adjust_trust_lvl_1`: (`float`) Adjusted Trust level 1 for next iteration.
            - `adjust_trust_lvl_2`: (`float`) Adjusted Trust level 2 for next iteration.
            
        """
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

    