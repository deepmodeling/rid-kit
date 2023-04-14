from typing import List, Optional, Dict
from pathlib import Path
import numpy as np
from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Parameter,
    BigParameter
)
from rid.utils import save_txt, set_directory
from rid.mcmc.walker import Walker,my_hist2d
from rid.select.model_devi import test_ef
from rid.common.tensorflow.graph import load_graph
try:
    import tensorflow.compat.v1 as tf
    tf.disable_v2_behavior()
except ImportError:
    import tensorflow as tf
from rid.constants import mcmc_cv_name
    
# kinetic enery in eV
kbT = (8.617343E-5) * 300 
beta = 1.0 / kbT

# ev to kj/mol
f_cvt = 96.485

class MCMCRun(OP):
    """
    `MCMC_Run` performs Markov chain Monte carlo sampling to reduce the dimension of Free energy surface (represented by neural network)
    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "task_names": str,
                "mcmc_config": BigParameter(Dict),
                "models": Artifact(Path, archive = None)
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "mcmc_cv": Artifact(Path, archive = None)
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

            - `task_name`: str,
            - `models`: Artifact(List[Path], optional=True),
          
        Returns
        -------
            Output dict with components:
        
            - `fes_fig`: (`Artifact(Path)`)
        """
        graph = load_graph(op_in["models"])
        mcmc_config = op_in["mcmc_config"]
        fd = mcmc_config["cv_dimension"]
        ns = mcmc_config["numb_steps"]
        nw = mcmc_config["numb_walkers"]
        cv1 = mcmc_config["cv_index_1"]
        cv2 = mcmc_config["cv_index_2"]
        cv_type = mcmc_config["cv_type"]
        bins = mcmc_config["bins"]
        
        if cv_type == "dih":
            xx = np.linspace(0,2* np.pi, bins)
            yy = np.linspace(0,2* np.pi, bins)
            pp_hist2d = np.zeros((1, len(xx), len(yy)))
            delta = 2.0 * np.pi / (bins-1)
        elif cv_type == "dis":
            xx = np.linspace(0,10, bins)
            yy = np.linspace(0,10, bins)
            pp_hist2d = np.zeros((1, len(xx), len(yy)))
            delta = 10.0 / (bins-1)
        else:
            raise ValueError("Undefined cv type, only support 'dih' and 'dis' type")
        
        task_path = Path(op_in["task_names"])
        task_path.mkdir(exist_ok=True, parents=True)

        with set_directory(task_path):
            with tf.Session(graph = graph) as sess:        
                walker = Walker(fd, nw, sess, cv_type)
                for _ in range(100):
                    pp, ee, ff = walker.sample(test_ef)

                for ii in range(ns+1):
                    pp, ee, ff = walker.sample(test_ef)
                    
                     ##certain 2d
                    pp_hist_new2d = my_hist2d(pp, xx, yy, delta, cv1, cv2)
                    pp_hist2d = (pp_hist2d * ii + pp_hist_new2d) / (ii+1)
                    
                    if ii == ns:
                        zz2d = np.transpose(-np.log(pp_hist2d+1e-10), (0,2,1))/beta
                        # convert ev to kcal/mol
                        zz2d *= f_cvt/4.184
                        zz2d = zz2d - np.min(zz2d)
                        np.savetxt(mcmc_cv_name,zz2d[0])
            
            
        op_out = OPIO(
            {
               "mcmc_cv": task_path.joinpath(mcmc_cv_name)
            }
        )
        return op_out