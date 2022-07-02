from rid.task.task import Task, TaskGroup
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Union
from rid.constants import (
        gmx_conf_name,
        gmx_top_name,
        gmx_mdp_name, 
        plumed_input_name,
        plumed_output_name,
        tf_graph_name
    )
from rid.utils import read_binary, write_binary, read_txt, write_txt
from rid.common.gromacs import make_md_mdp_from_config
from rid.common.plumed import make_deepfe_plumed, check_deepfe_input, make_restraint_plumed


class TaskBuilder(ABC):

    @abstractmethod
    def build(self):
        pass


class EnhcMDTaskBuilder(TaskBuilder):
    def __init__(
        self,
        conf: Optional[str],
        topology: Optional[str],
        gmx_config: Dict,
        cv_file: Optional[str] = None,
        selected_resid: Optional[List[int]] = None,
        trust_lvl_1: float = 1.0,
        trust_lvl_2: float = 2.0,
        model_list: List[str] = ["graph.pb"],
        plumed_output: str = "plm.out",
        cv_mode: str = "torsion"
    ):
        super().__init__()
        self.conf = conf
        self.topology = topology
        self.gmx_config = gmx_config
        self.stride = self.gmx_config["output_freq"]
        self.cv_file = cv_file
        self.selected_resid = selected_resid
        self.trust_lvl_1 = trust_lvl_1
        self.trust_lvl_2 = trust_lvl_2
        assert self.trust_lvl_1 < self.trust_lvl_2
        self.model_list = model_list
        self.plumed_output = plumed_output
        self.cv_mode = cv_mode
        self.task = Task()
    
    def build(self) -> Task:
        task_dict = {}
        task_dict.update(self.build_gmx())
        task_dict.update(self.build_plumed())
        # for idx, graph in enumerate(self.model_list):
        #     task_dict[tf_graph_name.format(idx=idx)] = (read_binary(graph), "wb")
        for fname, fconts in task_dict.items():
            self.task.add_file(fname, fconts)
        self.task.add_property({"num_models": len(self.model_list)})
        return self.task
    
    def get_task(self):
        return self.task
     
    def build_gmx(self):
        return build_gmx_dict(self.conf, self.topology, self.gmx_config)
        
    def build_plumed(self):
        return build_plumed_dict(
            conf=self.conf, cv_file=self.cv_file, selected_resid=self.selected_resid,
            trust_lvl_1=self.trust_lvl_1, trust_lvl_2=self.trust_lvl_2,
            model_list=self.model_list, stride=self.stride, output=self.plumed_output,
            mode=self.cv_mode
        )


class RestrainedMDTaskBuilder(TaskBuilder):
    def __init__(
        self,
        conf: Optional[str],
        topology: Optional[str],
        gmx_config: Dict,
        cv_file: Optional[str] = None,
        selected_resid: Optional[List[int]] = None,
        kappa: Union[int, float, List[Union[int, float]]] = 0.5,
        at: Union[int, float, List[Union[int, float]]] = 1.0,
        plumed_output: str = "plm.out",
        cv_mode: str = "torsion"
    ):
        super().__init__()
        self.conf = conf
        self.topology = topology
        self.gmx_config = gmx_config
        self.stride = self.gmx_config["output_freq"]
        self.cv_file = cv_file
        self.selected_resid = selected_resid
        self.plumed_output = plumed_output
        self.cv_mode = cv_mode
        self.kappa = kappa
        self.at = at
        self.task = Task()
    
    def build(self) -> Task:
        task_dict = {}
        task_dict.update(self.build_gmx())
        task_dict.update(self.build_plumed())
        for fname, fconts in task_dict.items():
            self.task.add_file(fname, fconts)
        return self.task
    
    def get_task(self):
        return self.task
     
    def build_gmx(self):
        return build_gmx_dict(self.conf, self.topology, self.gmx_config)
        
    def build_plumed(self):
        return build_plumed_restraint_dict(
            conf=self.conf, cv_file=self.cv_file, selected_resid=self.selected_resid,
            kappa=self.kappa, at=self.at,
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
    gmx_task_files[gmx_mdp_name]  = (make_md_mdp_from_config(gmx_config), "w")
    return gmx_task_files
    

def build_plumed_dict(
        conf: Optional[str] = None,
        cv_file: Optional[str] = None,
        selected_resid: Optional[List[int]] = None,
        trust_lvl_1: float = 1.0,
        trust_lvl_2: float = 2.0,
        model_list: List[str] = ["graph.pb"],
        stride: int = 100,
        output: str = "plm.out",
        mode: str = "torsion"    
    ):
    plumed_task_files = {}
    plm_content = make_deepfe_plumed(
        conf=conf, cv_file=cv_file, selected_resid=selected_resid,
        trust_lvl_1=trust_lvl_1, trust_lvl_2=trust_lvl_2,
        model_list=model_list, stride=stride,
        output=output, mode=mode
    )
    plumed_task_files[plumed_input_name] = (plm_content, "w")
    return plumed_task_files

def build_plumed_restraint_dict(
        conf: Optional[str] = None,
        cv_file: Optional[str] = None,
        selected_resid: Optional[List[int]] = None,
        kappa: Union[int, float, List[Union[int, float]]] = 0.5,
        at: Union[int, float, List[Union[int, float]]] = 1.0,
        stride: int = 100,
        output: str = "plm.out",
        mode: str = "torsion"    
    ):
    plumed_task_files = {}
    plm_content = make_restraint_plumed(
        conf=conf, cv_file=cv_file, selected_resid=selected_resid,
        kappa=kappa, at=at, stride=stride,
        output=output, mode=mode
    )
    plumed_task_files[plumed_input_name] = (plm_content, "w")
    return plumed_task_files