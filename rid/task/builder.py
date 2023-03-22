from rid.task.task import Task
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Union, Sequence
import numpy as np
from rid.constants import (
        gmx_conf_name,
        lmp_conf_name,
        gmx_top_name,
        lmp_input_name,
        gmx_mdp_name, 
        plumed_input_name,
    )
from rid.utils import read_txt
from rid.common.mol import get_distance_from_atomid
from rid.common.gromacs import make_md_mdp_string
from rid.common.plumed import make_deepfe_plumed, make_restraint_plumed, make_constraint_plumed, get_cv_name,make_distance_list_from_file


class TaskBuilder(ABC):

    @abstractmethod
    def build(self):
        pass


class EnhcMDTaskBuilder(TaskBuilder):
    def __init__(
        self,
        conf: str,
        topology: Optional[str],
        exploration_config: Dict,
        cv_file: Optional[List[str]] = None,
        selected_resid: Optional[List[int]] = None,
        selected_atomid: Optional[List[int]] = None,
        sampler_type: str = "gmx",
        trust_lvl_1: float = 1.0,
        trust_lvl_2: float = 2.0,
        model_list: List[str] = ["graph.pb"],
        plumed_output: str = "plm.out",
        cv_mode: str = "torsion",
        wall_list: Optional[List[str]] = None,
        iteration: Optional[str] = None
    ):
        super().__init__()
        self.conf = conf
        self.topology = topology
        self.exploration_config = exploration_config
        self.stride = self.exploration_config["output_freq"]
        self.cv_file = cv_file
        self.selected_resid = selected_resid
        self.selected_atomid = selected_atomid
        self.sampler_type = sampler_type
        self.trust_lvl_1 = trust_lvl_1
        self.trust_lvl_2 = trust_lvl_2
        assert self.trust_lvl_1 < self.trust_lvl_2
        self.model_list = model_list
        self.plumed_output = plumed_output
        self.cv_mode = cv_mode
        self.wall_list = wall_list
        self.iteration = iteration
        self.task = Task()
        self.cv_names = get_cv_name(
            conf=self.conf, cv_file=self.cv_file,
            selected_resid=self.selected_resid,
            selected_atomid=self.selected_atomid,
            stride=self.stride,
            mode=self.cv_mode
        )
    
    def build(self) -> Task:
        task_dict = {}
        if self.sampler_type == "gmx":
            task_dict.update(self.build_gmx())
        elif self.sampler_type == "lmp":
            task_dict.update(self.build_lmp())
        task_dict.update(self.build_plumed())
        for fname, fconts in task_dict.items():
            self.task.add_file(fname, fconts)
        self.task.add_property({"num_models": len(self.model_list)})
        return self.task
    
    def get_task(self):
        return self.task
     
    def build_gmx(self):
        return build_gmx_dict(self.conf, self.topology, self.exploration_config)
    
    def build_lmp(self):
        return build_lmp_dict(self.conf)
        
    def build_plumed(self):
        return build_plumed_dict(
            conf=self.conf, cv_file=self.cv_file, selected_resid=self.selected_resid,
            selected_atomid=self.selected_atomid,
            trust_lvl_1=self.trust_lvl_1, trust_lvl_2=self.trust_lvl_2,
            model_list=self.model_list, stride=self.stride, output=self.plumed_output,
            mode=self.cv_mode, wall_list = self.wall_list, iteration=self.iteration
        )
    
    def get_cv_dim(self):
        return len(self.cv_names)


class RestrainedMDTaskBuilder(TaskBuilder):
    def __init__(
        self,
        conf: str,
        topology: Optional[str],
        label_config: Dict,
        cv_file: Optional[List[str]] = None,
        selected_resid: Optional[List[int]] = None,
        selected_atomid: Optional[List[int]] = None,
        sampler_type: str = "gmx",
        kappa: Union[int, float, List[Union[int, float]]] = 0.5,
        at: Union[int, float, List[Union[int, float]]] = 1.0,
        plumed_output: str = "plm.out",
        cv_mode: str = "torsion"
    ):
        super().__init__()
        self.conf = conf
        self.topology = topology
        self.label_config = label_config
        self.stride = self.label_config["output_freq"]
        self.cv_file = cv_file
        self.selected_resid = selected_resid
        self.selected_atomid = selected_atomid
        self.plumed_output = plumed_output
        self.cv_mode = cv_mode
        self.sampler_type = sampler_type
        self.kappa = kappa
        self.at = at
        self.task = Task()
    
    def build(self) -> Task:
        task_dict = {}
        if self.sampler_type == "gmx":
            task_dict.update(self.build_gmx())
        elif self.sampler_type == "lmp":
            task_dict.update(self.build_lmp())
        task_dict.update(self.build_plumed())
        for fname, fconts in task_dict.items():
            self.task.add_file(fname, fconts)
        return self.task
    
    def get_task(self):
        return self.task
     
    def build_gmx(self):
        return build_gmx_dict(self.conf, self.topology, self.label_config)
    
    def build_lmp(self):
        return build_lmp_dict(self.conf)
        
    def build_plumed(self):
        return build_plumed_restraint_dict(
            conf=self.conf, cv_file=self.cv_file, selected_resid=self.selected_resid,
            selected_atomid=self.selected_atomid, kappa=self.kappa, at=self.at,
            stride=self.stride, output=self.plumed_output, mode=self.cv_mode
        )

class ConstrainedMDTaskBuilder(TaskBuilder):
    def __init__(
        self,
        conf: str,
        topology: Optional[str],
        label_config: Dict,
        cv_file: Optional[List[str]] = None,
        selected_atomid: Optional[List[int]] = None,
        sampler_type: str = "gmx",
        plumed_output: str = "plm.out",
        cv_mode: str = "distance"
    ):
        super().__init__()
        self.conf = conf
        self.topology = topology
        self.label_config = label_config
        self.stride = self.label_config["output_freq"]
        self.cv_file = cv_file
        self.selected_atomid = selected_atomid
        self.plumed_output = plumed_output
        self.cv_mode = cv_mode
        self.sampler_type = sampler_type
        self.task = Task()
    
    def build(self) -> Task:
        task_dict = {}
        if self.sampler_type == "gmx":
            task_dict.update(self.build_gmx())
        elif self.sampler_type == "lmp":
            task_dict.update(self.build_lmp())
        task_dict.update(self.build_plumed())
        for fname, fconts in task_dict.items():
            self.task.add_file(fname, fconts)
        return self.task
    
    def get_task(self):
        return self.task
     
    def build_gmx(self):
        return build_gmx_constraint_dict(self.conf, self.topology, self.label_config, self.selected_atomid)
    
    def build_lmp(self):
        return build_lmp_dict(self.conf)
        
    def build_plumed(self):
        return build_plumed_constraint_dict(
            conf=self.conf, cv_file=self.cv_file, selected_atomid=self.selected_atomid,
            stride=self.stride, output=self.plumed_output, mode=self.cv_mode
        )

def build_gmx_dict(
        conf: str,
        topology: str,
        gmx_config: Dict
    ):
    gmx_task_files = {}
    gmx_task_files[gmx_conf_name] = (read_txt(conf), "w")
    gmx_task_files[gmx_top_name]  = (read_txt(topology), "w")
    mdp_string = make_md_mdp_string(gmx_config)
    gmx_task_files[gmx_mdp_name]  = (mdp_string, "w")
    return gmx_task_files

def build_gmx_constraint_dict(
        conf: str,
        topology: str,
        gmx_config: Dict,
        selected_atomid: Optional[List[int]] = None
    ):
    gmx_task_files = {}
    gmx_task_files[gmx_conf_name] = (read_txt(conf), "w")
    cv_info = get_distance_from_atomid(conf, selected_atomid)
    ret = ""
    with open(topology, "r") as f:
        for line in f.readlines():
            ret += line
            if "constraints" in line:
                print("constrained md operating normally!\n")
                ret += "; atom1 atom2    funct   dis\n"
                for dis_id in range(len(selected_atomid)):
                    ret += "%s %s 2 %s\n"%(selected_atomid[dis_id][0], selected_atomid[dis_id][1],\
                        cv_info["%s %s"%(selected_atomid[dis_id][0],selected_atomid[dis_id][1])])        
    gmx_task_files[gmx_top_name]  = (ret, "w")
    mdp_string = make_md_mdp_string(gmx_config)
    gmx_task_files[gmx_mdp_name]  = (mdp_string, "w")
    return gmx_task_files

def build_lmp_dict(
    conf: str
):
    lmp_task_files = {}
    lmp_task_files[lmp_conf_name] = (read_txt(conf), "w")
    return lmp_task_files
    
def build_plumed_dict(
        conf: Optional[str] = None,
        cv_file: Optional[str] = None,
        selected_resid: Optional[List[int]] = None,
        selected_atomid: Optional[List[int]] = None,
        trust_lvl_1: float = 1.0,
        trust_lvl_2: float = 2.0,
        model_list: List[str] = ["graph.pb"],
        stride: int = 100,
        output: str = "plm.out",
        mode: str = "torsion",
        wall_list: Optional[List[str]] = None,
        iteration: Optional[str] = None
    ):
    plumed_task_files = {}
    plm_content = make_deepfe_plumed(
        conf=conf, cv_file=cv_file, selected_resid=selected_resid,
        selected_atomid = selected_atomid,
        trust_lvl_1=trust_lvl_1, trust_lvl_2=trust_lvl_2,
        model_list=model_list, stride=stride,
        output=output, mode=mode, wall_list=wall_list, iteration=iteration
    )
    plumed_task_files[plumed_input_name] = (plm_content, "w")
    return plumed_task_files

def build_plumed_restraint_dict(
        conf: Optional[str] = None,
        cv_file: Optional[str] = None,
        selected_resid: Optional[List[int]] = None,
        selected_atomid: Optional[List[int]] = None,
        kappa: Union[int, float, Sequence, np.ndarray] = 0.5,
        at: Union[int, float, Sequence, np.ndarray] = 1.0,
        stride: int = 100,
        output: str = "plm.out",
        mode: str = "torsion"
    ):
    plumed_task_files = {}
    if selected_atomid is not None:
        at = []
        cv_info = get_distance_from_atomid(conf, selected_atomid)
        for dis_id in range(len(selected_atomid)):
            at.append(cv_info["%s %s"%(selected_atomid[dis_id][0],selected_atomid[dis_id][1])])
    plm_content = make_restraint_plumed(
        conf=conf, cv_file=cv_file, selected_resid=selected_resid,selected_atomid = selected_atomid,
        kappa=kappa, at=at, stride=stride,
        output=output, mode=mode
    )
    plumed_task_files[plumed_input_name] = (plm_content, "w")
    return plumed_task_files

def build_plumed_constraint_dict(
        conf: Optional[str] = None,
        cv_file: Optional[str] = None,
        selected_atomid: Optional[List[int]] = None,
        stride: int = 100,
        output: str = "plm.out",
        mode: str = "distance"
    ):
    plumed_task_files = {}
    plm_content = make_constraint_plumed(
        conf=conf, cv_file=cv_file, selected_atomid=selected_atomid,
        stride=stride, output=output, mode=mode
    )
    plumed_task_files[plumed_input_name] = (plm_content, "w")
    return plumed_task_files