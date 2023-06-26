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
from rid.mcmc.walker import Walker, my_hist1d, my_hist2d, my_hist1d_path, my_hist2d_path
from rid.select.model_devi import test_ef
from rid.common.tensorflow.graph import load_graph
try:
    import tensorflow.compat.v1 as tf
    tf.disable_v2_behavior()
except ImportError:
    import tensorflow as tf
from rid.constants import mcmc_1cv_name, mcmc_1cv_dir_name, mcmc_2cv_name, kb, f_cvt, kcal2kj
import os

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
                "mcmc_1cv": Artifact(Path, archive = None),
                "mcmc_2cv": Artifact(Path, archive = None)
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
        temperature = mcmc_config["temperature"]
        fd = mcmc_config["cv_dimension"]
        ns = mcmc_config["numb_steps"]
        nw = mcmc_config["numb_walkers"]
        # kinetic enery in eV
        kbT = kb * temperature
        beta = 1.0 / kbT
        if "cv_upper_bound" in mcmc_config:
            cv_upper = mcmc_config["cv_upper_bound"]
        else:
            cv_upper = 2*np.pi
        if "cv_lower_bound" in mcmc_config:
            cv_lower = mcmc_config["cv_lower_bound"]
        else:
            cv_lower = 0
        proj_info = mcmc_config["proj_info"]
        proj_mode = proj_info["proj_mode"]
        proj_cv_index = proj_info["proj_cv_index"]
        if "path_list" in proj_info:
            path_list = np.array(proj_info["path_list"])
            path_lm = proj_info["path_lm"]
        cv_type = mcmc_config["cv_type"]
        bins = mcmc_config["bins"]
        
        if cv_type == "dih":
            xx = np.linspace(0,2* np.pi, bins)
            yy = np.linspace(0,2* np.pi, bins)
            pp_hist = np.zeros((fd, len(xx)))
            pp_hist2d = np.zeros((1, len(xx), len(yy)))
            delta = 2.0 * np.pi / (bins-1)
        elif cv_type == "dis":
            xx = np.linspace(0,10, bins)
            yy = np.linspace(0,10, bins)
            if proj_mode == "cv":
                pp_hist = np.zeros((fd, len(xx)))
            elif proj_mode == "path":
                pp_hist = np.zeros((2, len(xx)))
            pp_hist2d = np.zeros((1, len(xx), len(yy)))
            delta = 10.0 / (bins-1)
        else:
            raise ValueError("Undefined cv type, only support 'dih' and 'dis' type")
        
        task_path = Path(op_in["task_names"])
        task_path.mkdir(exist_ok=True, parents=True)

        mcmc_2cv_path = None
        with set_directory(task_path):
            with tf.Session(graph = graph) as sess:        
                walker = Walker(fd, nw, sess, cv_type, cv_lower=cv_lower, cv_upper=cv_upper)
                for _ in range(10000):
                    pp, ee, ff = walker.sample(test_ef)

                for ii in range(ns+1):
                    pp, ee, ff = walker.sample(test_ef)
                        
                    if proj_mode == "cv":
                        # project on 1D CV
                        pp_hist_new = my_hist1d(pp, xx, delta, fd)
                        pp_hist = (pp_hist * ii + pp_hist_new) / (ii+1)
                        if not os.path.exists(mcmc_1cv_dir_name):
                            os.makedirs(mcmc_1cv_dir_name)
                        if np.mod(ii,int(ns/5)) == 0:
                            zz = -np.log(pp_hist+1e-7)/beta
                            # convert ev to kcal/mol
                            zz *= f_cvt/kcal2kj
                            zz = zz - np.min(zz)      
                            for jj in range(fd):
                                fp = open(mcmc_1cv_dir_name+"/"+mcmc_1cv_name.format(tag=jj), "a")
                                for temp in zz[jj]:
                                    fp.write(str(temp)+'    ')
                                fp.write('\n')
                                fp.close()
                        # project on 2D CV
                        assert len(proj_cv_index) == 2
                        cv1 = proj_cv_index[0]
                        cv2 = proj_cv_index[1]
                        ##certain 2d
                        pp_hist_new2d = my_hist2d(pp, xx, yy, delta, cv1, cv2)
                        pp_hist2d = (pp_hist2d * ii + pp_hist_new2d) / (ii+1)
                        if ii == ns:
                            zz2d = np.transpose(-np.log(pp_hist2d+1e-10), (0,2,1))/beta
                            # convert ev to kcal/mol
                            zz2d *= f_cvt/kcal2kj
                            zz2d = zz2d - np.min(zz2d)
                            np.savetxt(mcmc_2cv_name,zz2d[0])
                    elif proj_mode == "path":
                        # project on 1D CV
                        pp_hist_new = my_hist1d_path(pp, xx, delta, path_lm, path_list, proj_cv_index)
                        pp_hist = (pp_hist * ii + pp_hist_new) / (ii+1)
                        if not os.path.exists(mcmc_1cv_dir_name):
                            os.makedirs(mcmc_1cv_dir_name)
                        if np.mod(ii,int(ns/5)) == 0:
                            zz = -np.log(pp_hist+1e-7)/beta
                            # convert ev to kcal/mol
                            zz *= f_cvt/kcal2kj
                            zz = zz - np.min(zz)    
                            # iterate over 2 path CV  
                            for jj in range(2):
                                fp = open(mcmc_1cv_dir_name+"/"+mcmc_1cv_name.format(tag=jj), "a")
                                for temp in zz[jj]:
                                    fp.write(str(temp)+'    ')
                                fp.write('\n')
                                fp.close()
                        # project on 2D CV
                        pp_hist_new2d_path = my_hist2d_path(pp, xx, yy, delta, path_lm, path_list, proj_cv_index)
                        pp_hist2d = (pp_hist2d * ii + pp_hist_new2d_path) / (ii+1)
                        if ii == ns:
                            zz2d = np.transpose(-np.log(pp_hist2d+1e-10), (0,2,1))/beta
                            # convert ev to kcal/mol
                            zz2d *= f_cvt/kcal2kj
                            zz2d = zz2d - np.min(zz2d)
                            np.savetxt(mcmc_2cv_name,zz2d[0])
                            
        if os.path.exists(task_path.joinpath(mcmc_2cv_name)):
            mcmc_2cv_path = task_path.joinpath(mcmc_2cv_name)
            
            
        op_out = OPIO(
            {
               "mcmc_1cv": task_path.joinpath(mcmc_1cv_dir_name),
               "mcmc_2cv": mcmc_2cv_path
            }
        )
        return op_out