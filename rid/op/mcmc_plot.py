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
from rid.constants import mcmc_1cv_dir_name, mcmc_1cv_name, mcmc_2cv_name,mcmc_1cv_fig, mcmc_2cv_fig, mcmc_2cv_fig_separate
from matplotlib import pyplot as plt
import os
    
# kinetic enery in eV
kbT = (8.617343E-5) * 300 
beta = 1.0 / kbT

# ev to kj/mol
f_cvt = 96.485

class MCMCPlot(OP):
    """
    `MCMC_Plot` plot the reduced free energy surface produced by MCMC run.
    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "mcmc_1cv": Artifact(List[Path]),
                "mcmc_2cv": Artifact(List[Path]),
                "mcmc_config": BigParameter(Dict)
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "mcmc_fig": Artifact(Path, archive = None)
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

            - `mcmc_config`: Dict,
            - `mcmc_cv`: Artifact(List[Path]),
          
        Returns
        -------
            Output dict with components:
        
            - `mcmc_fig`: (`Artifact(Path)`)
        """
        mcmc_config = op_in["mcmc_config"]
        cv_type = mcmc_config["cv_type"]
        bins = mcmc_config["bins"]
        cv_dim = mcmc_config["cv_dimension"]
        
        task_path = Path("mcmc_fig")
        task_path.mkdir(exist_ok=True, parents=True)
        with set_directory(task_path):
            fourdata_1cv = []
            for dir in op_in["mcmc_1cv"]:
                all_1cv = []
                for cv_index in range(cv_dim):
                    file = dir/mcmc_1cv_name.format(tag = cv_index)
                    pmf = np.loadtxt(file)
                    pmf = pmf - np.min(pmf)
                    proj_iterations = pmf.shape[0]
                    all_1cv.append(pmf)
                fourdata_1cv.append(all_1cv)
            avedata_1cv = np.mean(np.array(fourdata_1cv),axis=0)
            print("avedata shape", avedata_1cv.shape)
            # make 1cv plot
            if not os.path.exists(mcmc_1cv_dir_name):
                os.makedirs(mcmc_1cv_dir_name)
            if cv_type == "dih":
                xedges = np.linspace(0, 2*np.pi, bins)
            elif cv_type == "dis":
                xedges = np.linspace(0, 10, bins)
            for cv_index in range(cv_dim):
                plt.figure(figsize=(8, 6))
                for proj_iter in range(proj_iterations):
                    print("xedges shape", xedges.shape)
                    plt.plot(xedges,avedata_1cv[cv_index][proj_iter], label = proj_iter)
                plt.xlabel(r'cv')
                plt.ylabel(r'free energy (kcal/mol)')
                plt.legend()
                plt.savefig(mcmc_1cv_dir_name+"/"+mcmc_1cv_fig.format(tag = cv_index),dpi=600,bbox_inches='tight')
                
            fourdata_2cv = []
            for file in op_in["mcmc_2cv"]:
                pmf = np.loadtxt(file)
                pmf = pmf - np.min(pmf)
                allda=[]
                if cv_type == "dih":
                    xedges = np.linspace(0, 2*np.pi, bins)
                    yedges = np.linspace(0, 2*np.pi, bins)
                elif cv_type == "dis":
                    xedges = np.linspace(0, 10, bins)
                    yedges = np.linspace(0, 10, bins)
                for i in range(bins):
                    for j in range(bins):
                        temp = []
                        temp.append(xedges[i])
                        temp.append(yedges[j])
                        temp.append(pmf[i][j])
                        temp = np.array(temp)
                        allda.append(temp)
                newarray=np.array(allda)
                idex=np.lexsort([newarray[:,1], newarray[:,0]])
                sorted_data = newarray[idex,:]
                fourdata_2cv.append(sorted_data[:,2])
            avedata_2cv = np.mean(np.array(fourdata_2cv),axis=0)    
            # make 2cv plot
            if cv_type == "dih":
                xedges = np.linspace(0, 2*np.pi, bins)
                yedges = np.linspace(0, 2*np.pi, bins)
            elif cv_type == "dis":
                xedges = np.linspace(0, 10, bins)
                yedges = np.linspace(0, 10, bins)
            plt.figure(figsize=(8, 6))
            cmap = plt.cm.get_cmap("jet_r")
            CS = plt.contourf(xedges,yedges,avedata_2cv.reshape(bins,bins),levels = np.linspace(0,6,61),cmap=cmap,extend="max")
            cbar = plt.colorbar(CS)
            cbar.ax.tick_params(labelsize=8) 
            cbar.ax.set_title('kcal/mol',fontsize=8)
            plt.xlabel(r'CV index 1')
            plt.ylabel(r'CV index 2')
            plt.savefig(mcmc_2cv_fig,dpi=600,bbox_inches='tight')
            for iii in range(len(fourdata_2cv)):
                fig = plt.figure(figsize=(8, 6))
                cmap = plt.cm.get_cmap("jet_r")
                CS = plt.contourf(xedges,yedges,fourdata_2cv[iii].reshape(bins,bins),levels = np.linspace(0,6,61),cmap=cmap,extend="max")
                cbar = plt.colorbar(CS)
                cbar.ax.tick_params(labelsize=8) 
                cbar.ax.set_title('kcal/mol',fontsize=8)
                plt.xlabel(r'CV index 1')
                plt.ylabel(r'CV index 2')
                plt.savefig(mcmc_2cv_fig_separate.format(tag=iii),dpi=600,bbox_inches='tight')
                
        op_out = OPIO(
            {
               "mcmc_fig": task_path
            }
        )
        return op_out