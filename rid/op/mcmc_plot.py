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
from rid.constants import mcmc_cv_name, mcmc_cv_fig, mcmc_cv_fig_separate
from matplotlib import pyplot as plt
    
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
                "mcmc_cv": Artifact(List[Path]),
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
        
        task_path = Path("mcmc_fig")
        task_path.mkdir(exist_ok=True, parents=True)
        with set_directory(task_path):
            fourdata=[]
            for file in op_in["mcmc_cv"]:
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
                fourdata.append(sorted_data[:,2])
            avedata=np.mean(np.array(fourdata),axis=0)
            if cv_type == "dih":
                xedges = np.linspace(0, 2*np.pi, bins)
                yedges = np.linspace(0, 2*np.pi, bins)
            elif cv_type == "dis":
                xedges = np.linspace(0, 10, bins)
                yedges = np.linspace(0, 10, bins)
            fig = plt.figure(figsize=(8, 6))
            cmap = plt.cm.get_cmap("jet_r")
            CS = plt.contourf(xedges,yedges,avedata.reshape(bins,bins),levels = np.linspace(0,6,61),cmap=cmap,extend="max")
            cbar = plt.colorbar(CS)
            plt.xlabel(r'a')
            plt.ylabel(r'b')
            plt.savefig(mcmc_cv_fig,dpi=600,bbox_inches='tight')
            for iii in range(len(fourdata)):
                fig = plt.figure(figsize=(8, 6))
                cmap = plt.cm.get_cmap("jet_r")
                CS = plt.contourf(xedges,yedges,fourdata[iii].reshape(bins,bins),levels = np.linspace(0,6,61),cmap=cmap,extend="max")
                cbar = plt.colorbar(CS)
                plt.xlabel(r'a')
                plt.ylabel(r'b')
                plt.savefig(mcmc_cv_fig_separate.format(tag=iii),dpi=600,bbox_inches='tight')
                
        op_out = OPIO(
            {
               "mcmc_fig": task_path
            }
        )
        return op_out