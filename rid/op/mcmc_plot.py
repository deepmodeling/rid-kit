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
                "plm_out": Artifact(List[Path], optional=True),
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
        proj_info = mcmc_config["proj_info"]
        proj_mode = proj_info["proj_mode"]
        
        if proj_mode == "cv":
            proj_1d_num = cv_dim
        elif proj_mode == "path":
            proj_1d_num = 2
        else:
            raise ValueError("Invalid cv type")
        
        task_path = Path("mcmc_fig")
        task_path.mkdir(exist_ok=True, parents=True)
        with set_directory(task_path):
            fourdata_1cv = []
            for dir in op_in["mcmc_1cv"]:
                all_1cv = []
                for cv_index in range(proj_1d_num):
                    file = dir/mcmc_1cv_name.format(tag = cv_index)
                    pmf = np.loadtxt(file)
                    pmf = pmf - np.min(pmf)
                    allda=[]
                    proj_iterations = pmf.shape[0]
                    if cv_type == "dih":
                        xedges = np.linspace(0, 2*np.pi, bins)
                        for i in range(bins):
                            temp = []
                            # this is to change to dihedral dimension to (-pi,pi) which corresponds to the gmx output
                            if xedges[i]>=np.pi:
                                temp.append(xedges[i] - np.pi*2)
                            else:
                                temp.append(xedges[i] + np.pi*2/(bins-1))
                            if i == (bins - 1):
                                temp.append(pmf[:,0])
                            else:
                                temp.append(pmf[:,i])
                            temp = np.array(temp)
                            allda.append(temp)
                    elif cv_type == "dis":
                        xedges = np.linspace(0, 10, bins)
                        for i in range(bins):
                            temp = []
                            temp.append(xedges[i])
                            temp.append(pmf[:,i])
                            temp = np.array(temp)
                            allda.append(temp)
                    newarray=np.array(allda)
                    idex=np.argsort(newarray[:,0])
                    sorted_data = newarray[idex,:]
                    print("sorted data shape",sorted_data.shape)
                    sorted_data = np.stack(sorted_data[:,1])
                    all_1cv.append(sorted_data)
                all_1cv = np.array(all_1cv)
                print("all_1cv shape", all_1cv.shape)
                fourdata_1cv.append(all_1cv)
            avedata_1cv = np.mean(np.array(fourdata_1cv),axis=0)
            print("avedata shape", avedata_1cv.shape)
            # make 1cv plot
            if not os.path.exists(mcmc_1cv_dir_name):
                os.makedirs(mcmc_1cv_dir_name)

            for cv_index in range(proj_1d_num):
                if proj_mode == "cv":
                    if cv_type == "dih":
                        xedges = np.linspace(-np.pi, np.pi, bins)
                    elif cv_type == "dis":
                        xedges = np.linspace(0, 10, bins)
                    plt.figure(figsize=(8, 6))
                    for proj_iter in range(proj_iterations):
                        plt.plot(xedges,avedata_1cv[cv_index][:,proj_iter], label = proj_iter)
                    plt.xlabel(r'cv')
                    plt.ylabel(r'free energy (kcal/mol)')
                    plt.legend()
                    plt.savefig(mcmc_1cv_dir_name+"/"+mcmc_1cv_fig.format(tag = cv_index),dpi=600,bbox_inches='tight')
                elif proj_mode == "path":
                    path_edges = [np.linspace(0, 10, bins),np.linspace(-10, 0, bins)]
                    plt.figure(figsize=(8, 6))
                    for proj_iter in range(proj_iterations):
                        plt.plot(path_edges[cv_index],avedata_1cv[cv_index][:,proj_iter], label = proj_iter)
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
                    for i in range(bins):
                        for j in range(bins):
                            temp = []
                            # this is to change to dihedral dimension to (-pi,pi) which corresponds to the gmx output
                            if xedges[i]>=np.pi:
                                temp.append(xedges[i] - np.pi*2)
                            else:
                                temp.append(xedges[i] + np.pi*2/(bins-1))
                            if yedges[j]>=np.pi:
                                temp.append(yedges[j] - np.pi*2)
                            else:
                                temp.append(yedges[j] + np.pi*2/(bins-1))
                            if j == (bins -1):
                                temp.append(pmf[i][0])
                            elif i == (bins - 1):
                                temp.append(pmf[0][j])
                            else:
                                temp.append(pmf[i][j])
                            temp = np.array(temp)
                            allda.append(temp)
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
                xedges = np.linspace(-np.pi, np.pi, bins)
                yedges = np.linspace(-np.pi, np.pi, bins)
            elif cv_type == "dis":
                if proj_mode == "cv":
                    xedges = np.linspace(0, 10, bins)
                    yedges = np.linspace(0, 10, bins)
                elif proj_mode == "path":
                    xedges = np.linspace(0, 10, bins)
                    yedges = np.linspace(-10, 0, bins)
            plt.figure(figsize=(8, 6))
            cmap = plt.cm.get_cmap("jet_r")
            # Define percentiles for the levels
            upper_perc = np.percentile(np.unique(fourdata_2cv[0]), 99)
            CS = plt.contourf(xedges,yedges,avedata_2cv.reshape(bins,bins),levels = np.linspace(0,upper_perc,101),cmap=cmap,extend="max")
            
            if op_in["plm_out"] is not None:
                for cv_output in op_in["plm_out"]:
                    cv_point = np.loadtxt(cv_output)
                    assert cv_point.shape[1] == 2
                    point_numbers = list(range(len(cv_point[1:,0])))
                    P1 = plt.scatter(cv_point[1:,0], cv_point[1:,1], s = 2, marker = 'o', c = point_numbers)
                    P_init = plt.scatter(cv_point[0,0], cv_point[0,1], s = 20, marker = 'x', c = "k")
            cbar = plt.colorbar(CS)
            cbar.ax.tick_params(labelsize=8) 
            cbar.ax.set_title('kcal/mol',fontsize=8)
            if op_in["plm_out"] is not None:
                plt.legend([P_init],["init"])
                cbar2 = plt.colorbar(P1)
                cbar2.ax.tick_params(labelsize=4) 
                cbar2.ax.set_title('steps',fontsize=6)
            plt.xlabel(r'CV index 1')
            plt.ylabel(r'CV index 2')
            plt.savefig(mcmc_2cv_fig,dpi=600,bbox_inches='tight')
            for iii in range(len(fourdata_2cv)):
                fig = plt.figure(figsize=(8, 6))
                cmap = plt.cm.get_cmap("jet_r")
                # Define percentiles for the levels
                upper_perc = np.percentile(np.unique(fourdata_2cv[iii]), 99)
                CS = plt.contourf(xedges,yedges,fourdata_2cv[iii].reshape(bins,bins),levels = np.linspace(0,upper_perc,101),cmap=cmap,extend="max")
                if op_in["plm_out"] is not None:
                    for cv_output in op_in["plm_out"]:
                        cv_point = np.loadtxt(cv_output)
                        assert cv_point.shape[1] == 2
                        point_numbers = list(range(len(cv_point[1:,0])))
                        P1 = plt.scatter(cv_point[1:,0], cv_point[1:,1], s = 2, marker = 'o', c = point_numbers)
                        P_init = plt.scatter(cv_point[0,0], cv_point[0,1], s = 20, marker = 'x', c = "k")
                cbar = plt.colorbar(CS)
                cbar.ax.tick_params(labelsize=8) 
                cbar.ax.set_title('kcal/mol',fontsize=8)
                if op_in["plm_out"] is not None:
                    plt.legend([P_init],["init"])
                    cbar2 = plt.colorbar(P1)
                    cbar2.ax.tick_params(labelsize=4) 
                    cbar2.ax.set_title('steps',fontsize=6)
                plt.xlabel(r'CV index 1')
                plt.ylabel(r'CV index 2')
                plt.savefig(mcmc_2cv_fig_separate.format(tag=iii),dpi=600,bbox_inches='tight')
                
        op_out = OPIO(
            {
               "mcmc_fig": task_path
            }
        )
        return op_out