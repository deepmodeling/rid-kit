from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Parameter
)

from typing import Tuple, List, Optional, Dict
from pathlib import Path
from rid.constants import (
        explore_task_pattern, 
        gmx_conf_name,
        gmx_top_name,
        gmx_mdp_name, 
        plumed_input_name,
        plumed_output_name
    )
from rid.task.builder import EnhcMDTaskBuilder


class PrepExplore(OP):

    r"""Prepare files for exploration tasks.
    Currently, RiD is based on Gromacs with PLUMED2 plugin. Input files must contain .gro for conformations,
    .top for topology (with pre-defined force fields) and corresponding configuration in Dict formats while `models` 
    (graph files, .pb) are optional. Exploration step would run biased MD sampling with neural network models or 
    brute force MD sampling without neural network model provided.

    With models provided, the bias forces will be the average value of outputs of these models and tuned by a switching function.
    .. math::
        F(r) = -\nabla_{r_i} U(r) + \sigma( \me ( s( r))) \nabla_{r_i} A(r)
    where :math:`F(r)` is forces exrted on atoms, :math:`U(r)` is potential energy and :math:`A(r)` is free energy 
    represented by neural networks.
    .. math::
        \sigma(\epsilon)=
            \begin{cases}
                    1, & \epsilon<\epsilon_0 \\
                    \frac{1}{2}+\frac{1}{2}\cos{(\pi \frac{\epsilon-\epsilon_0}{\epsilon_1-\epsilon_0})}, & \epsilon_0 <\epsilon < \epsilon_1 \\
                    0, &\epsilon > \epsilon_1
            \end{cases}
    where :math:`\sigma(\epsilon)` is the switching function with parameters trust level (`trust_lvl_1` and `trust_lvl_2`).

    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "models": Artifact(List[Path], optional=True),
                "topology": Artifact(Path),
                "conf": Artifact(Path),
                "trust_lvl_1": float,
                "trust_lvl_2": float,
                "gmx_config": Dict,
                "cv_config": Dict,
                "task_name": str,
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "task_path": Artifact(Path),
                "cv_dim": int
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
        
            - `models`: (`Artifact(List[Path])`) Optional. Neural network model files (`.pb`) used to bias the simulation. 
                Run brute force MD simulations if not provided.
            - `trust_lvl_1`: (`float`) Trust level 1.
            - `trust_lvl_2`: (`float`) Trust level 2.
            - `topology`: (`Artifact(Path)`) Topology files (.top) for Gromacs simulations.
            - `conf`: (`Artifact(Path)`) Conformation files (.gro) for Gromacs simulations.
            - `gmx_config`: (`Dict`) Configuration in `Dict` format for Gromacs run. Must contains:
                `dt`, `steps`, `temperature`, `output_freq`.
            - `cv_config`: (`Dict`) Configuration for CV creation.
            - `task_name`: (`str`) Task name used to make sub-dir for tasks.
           
        Returns
        -------
            Output dict with components:
        
            - `task_path`: (`Artifact(Path)`) A directory containing files for RiD exploration.
            - `cv_dim`: (`int`) CV dimensions.
        """

        if op_in["cv_config"]["mode"] == "torsion":
            cv_file = None,
            selected_resid = op_in["cv_config"]["selected_resid"]
        elif op_in["cv_config"]["mode"] == "custom":
            cv_file = op_in["cv_config"]["cv_file"],
            selected_resid = None
        if op_in["models"] is None:
            models = []
        else:
            models = [str(model) for model in op_in["models"]]

        gmx_task_builder = EnhcMDTaskBuilder(
            conf = op_in["conf"],
            topology = op_in["topology"],
            gmx_config = op_in["gmx_config"],
            cv_file = cv_file,
            selected_resid = selected_resid,
            trust_lvl_1 = op_in["trust_lvl_1"],
            trust_lvl_2 = op_in["trust_lvl_2"],
            model_list = models,
            plumed_output = plumed_output_name,
            cv_mode = op_in["cv_config"]["mode"]
        )
        gmx_task = gmx_task_builder.build()
        cv_dim = gmx_task_builder.get_cv_dim()
        task_path = Path(op_in["task_name"])
        task_path.mkdir(exist_ok=True, parents=True)
        for fname, fconts in gmx_task.files.items():
            with open(task_path.joinpath(fname), fconts[1]) as ff:
                ff.write(fconts[0])
        op_out = OPIO(
            {
                "task_path": task_path,
                "cv_dim": int(cv_dim)
            }
        )
        return op_out