from logging import raiseExceptions
from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    BigParameter
)

from typing import List, Dict, Union
from pathlib import Path
from rid.constants import (
        plumed_output_name
    )
from rid.task.builder import EnhcMDTaskBuilder


class PrepExplore(OP):

    r"""
    Prepare files for exploration tasks.
    Currently, RiD is based on Gromacs/Lammps with PLUMED2 plugin. Provide .gro files and .top files if running Gromacs and .lmp files if running Lammps.
    Exploration step would run biased MD sampling with neural network models or 
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
                "topology": Artifact(Path, optional=True),
                "conf": Artifact(Path),
                "cv_file": Artifact(List[Path], optional=True),
                "trust_lvl_1": float,
                "trust_lvl_2": float,
                "exploration_config": BigParameter(Dict),
                "cv_config": BigParameter(Dict),
                "task_name": str,
                "block_tag": str
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "task_path": Artifact(Path, archive = None),
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
            - `conf`: (`Artifact(Path)`) Conformation files (.gro, .lmp) for Gromacs/Lammps simulations.
            - `exploration_config`: (`Dict`) Configuration in `Dict` format for Gromacs/Lammps run.
            - `cv_config`: (`Dict`) Configuration for CV creation.
            - `task_name`: (`str`) Task name used to make sub-dir for tasks.
           
        Returns
        -------
            Output dict with components:
        
            - `task_path`: (`Artifact(Path)`) A directory containing files for RiD exploration.
            - `cv_dim`: (`int`) CV dimensions.
        """
        cv_file = []
        selected_resid = None
        selected_atomid = None
        if op_in["cv_config"]["mode"] == "torsion":
            selected_resid = op_in["cv_config"]["selected_resid"]
        elif op_in["cv_config"]["mode"] == "distance":
            selected_atomid = op_in["cv_config"]["selected_atomid"]
        elif op_in["cv_config"]["mode"] == "custom":
            cv_file = op_in["cv_file"]
        if op_in["models"] is None:
            models = []
        else:
            models = [str(model.name) for model in op_in["models"]]

        wall_list = None
        if "iterative_walls" in op_in["cv_config"]:
            wall_list = op_in["cv_config"]["iterative_walls"]
        
        iteration = int(op_in["block_tag"].split("-")[1])
        gmx_task_builder = EnhcMDTaskBuilder(
            conf = op_in["conf"],
            topology = op_in["topology"],
            exploration_config = op_in["exploration_config"],
            cv_file=cv_file,
            selected_resid = selected_resid,
            selected_atomid = selected_atomid,
            sampler_type = op_in["exploration_config"]["type"],
            trust_lvl_1 = op_in["trust_lvl_1"],
            trust_lvl_2 = op_in["trust_lvl_2"],
            model_list = models,
            plumed_output = plumed_output_name,
            cv_mode = op_in["cv_config"]["mode"],
            wall_list = wall_list,
            iteration = iteration
        )
        cv_dim = gmx_task_builder.get_cv_dim()
        task_path = Path(op_in["task_name"])
        task_path.mkdir(exist_ok=True, parents=True)
        gmx_task = gmx_task_builder.build()
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