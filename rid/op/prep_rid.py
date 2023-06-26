import os, sys, shutil, logging
from typing import List, Dict
from pathlib import Path
from copy import deepcopy
from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    BigParameter
)
from rid.utils import load_json
from rid.constants import model_tag_fmt, init_conf_gmx_name, init_conf_lmp_name,init_input_name, walker_tag_fmt


logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def prep_confs(confs, numb_walkers,sampler_type):
    numb_confs = len(confs)
    conf_list = []
    if sampler_type == "gmx":
        if numb_confs < numb_walkers:
            logger.info("Number of confs is smaller than number of walkers. Copy replicas up to number of walkers.")
            for idx in range(numb_walkers):
                shutil.copyfile(confs[idx%numb_confs], init_conf_gmx_name.format(idx=idx))
        elif numb_confs > numb_walkers:
            logger.info("Number of confs is greater than number of walkers. Only use the fist `numb_walkers` confs.")
            for idx in range(numb_walkers):
                shutil.copyfile(confs[idx], init_conf_gmx_name.format(idx=idx))
        else:
            for idx in range(numb_walkers):
                shutil.copyfile(confs[idx], init_conf_gmx_name.format(idx=idx))
        for idx in range(numb_walkers):
            conf_list.append(Path(init_conf_gmx_name.format(idx=idx)))
    elif sampler_type == "lmp":
        if numb_confs < numb_walkers:
            logger.info("Number of confs is smaller than number of walkers. Copy replicas up to number of walkers.")
            for idx in range(numb_walkers):
                shutil.copyfile(confs[idx%numb_confs], init_conf_lmp_name.format(idx=idx))
        elif numb_confs > numb_walkers:
            logger.info("Number of confs is greater than number of walkers. Only use the fist `numb_walkers` confs.")
            for idx in range(numb_walkers):
                shutil.copyfile(confs[idx], init_conf_lmp_name.format(idx=idx))
        else:
            for idx in range(numb_walkers):
                shutil.copyfile(confs[idx], init_conf_lmp_name.format(idx=idx))
        for idx in range(numb_walkers):
            conf_list.append(Path(init_conf_lmp_name.format(idx=idx)))
    else:
        raise ValueError("Invalid sampler type, only support gmx and lmp")
    return conf_list


class PrepRiD(OP):

    """Pre-processing of RiD.
    
    1. Parse RiD configuration JSON file, get default value if parameters are not provided.
    2. Rearrange conformation files.
    3. Make task names and formats.
    """

    @classmethod
    def get_input_sign(cls):
        return OPIOSign(
            {
                "confs": Artifact(List[Path]),
                "rid_config": Artifact(Path)
            }
        )

    @classmethod
    def get_output_sign(cls):
        return OPIOSign(
            {
                "numb_iters": int,
                "numb_walkers": int,
                "numb_models": int,
                "confs": Artifact(List[Path],archive = None),
                "walker_tags": List,
                "model_tags": List,
                
                "exploration_config": BigParameter(Dict),
                "cv_config": BigParameter(Dict),
                "trust_lvl_1": List[float],
                "trust_lvl_2": List[float],
                "cluster_threshold": List[float],
                "angular_mask": List,
                "weights": List,
                "numb_cluster_upper": int,
                "numb_cluster_lower": int,
                "max_selection": int,
                "numb_cluster_threshold": int,
                "std_threshold": float,
                "dt": float,
                "output_freq": float,
                "slice_mode": str,
                "type_map": List,
                "label_config": BigParameter(Dict),
                "train_config": BigParameter(Dict)
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
        
            - `confs`: (`Artifact(Path)`) User-provided initial conformation files (.gro) for reinfoced dynamics.
            - `rid_config`: (`Artifact(Path)`) Configuration file (.json) of RiD. 
                Parameters in this file will be parsed.
           
        Returns
        -------
            Output dict with components:
        
            - `numb_iters`: (`int`) Max number of iterations for Ri.
            - `numb_walkers`: (`int`) Number of parallel walkers for exploration.
            - `numb_models`: (`int`) Number of neural network models of RiD.
            - `confs`: (`Artifact(List[Path])`) Rearranged initial conformation files (.gro) for reinfoced dynamics.
            - `walker_tags`: (`List`) Tag formmat for parallel walkers.
            - `model_tags`: (`List`) Tag formmat for neural network models.
            
            - `exploration_config`: (`Dict`) Configuration of simulations in exploration steps.
            - `cv_config`: (`Dict`) Configuration to create CV in PLUMED2 style.
            - `trust_lvl_1`: (`List[float]`) Trust level 1, or e0.
            - `trust_lvl_2`: (`List[float]`) Trust level 2, or e1.
            - `cluster_threshold`: (`List[float]`) Initial guess of cluster threshold.
            - `angular_mask`: (`List`) Angular mask for periodic collective variables. 
                1 represents periodic, 0 represents non-periodic.
            - `weights`: (`List`) Weights for clustering collective variables. see details in cluster algorithms.
            - `numb_cluster_upper`: (`int`) Upper limit of cluster number to make cluster threshold.
            - `numb_cluster_lower`: (`int`) Lower limit of cluster number to make cluster threshold.
            - `max_selection`: (`int`) Max selection number of clusters in Selection steps for each parallel walker.
            - `numb_cluster_threshold`: (`int`) Used to adjust trust level. When cluster number is grater than this threshold, 
                trust levels will be increased adaptively.
            - `dt`: (`float`) Time interval of exploration MD simulations. Gromacs `trjconv` commands will need this parameters 
                to slice trajectories by `-dump` tag, see `selection` steps for detail.
            - `slice_mode`: (`str`) Mode to slice trajectories. Either `gmx` or `mdtraj`.
            - `label_config`: (`Dict`) Configuration of simulations in labeling steps.
            - `kappas`: (`List`) Force constants of harmonic restraints to perform restrained MD simulations.
            - `train_config`: (`Dict`) Configuration to train neural networks, including training strategy and network structures.
        """

        jdata = deepcopy(load_json(op_in["rid_config"]))
        numb_walkers = jdata.pop("numb_walkers")
        train_config = jdata.pop("Train")
        numb_models = train_config.pop("numb_models")
        numb_iters = jdata.pop("numb_iters")
        exploration_config = jdata.pop("ExploreMDConfig")
        
        sampler_type = exploration_config["type"]
        conf_list = prep_confs(op_in["confs"], numb_walkers, sampler_type)

        walker_tags = []
        model_tags = []
        for idx in range(numb_walkers):
            walker_tags.append(walker_tag_fmt.format(idx=idx))
        for idx in range(numb_models):
            model_tags.append(model_tag_fmt.format(idx=idx))
        
        dt = exploration_config["dt"]
        output_freq = exploration_config["output_freq"]
        cv_config = jdata.pop("CV")
        angular_mask = cv_config["angular_mask"]
        weights = cv_config["weights"]
        
        selection_config = jdata.pop("SelectorConfig")

        label_config = jdata.pop("LabelMDConfig")
        
        trust_lvl_1 = jdata.pop("trust_lvl_1")
        trust_lvl_2 = jdata.pop("trust_lvl_2")
        trust_lvl_1_list = [trust_lvl_1 for _ in range(numb_walkers)]
        trust_lvl_2_list = [trust_lvl_2 for _ in range(numb_walkers)]

        cluster_threshold = selection_config.pop("cluster_threshold")
        cluster_threshold_list = [cluster_threshold for _ in range(numb_walkers)]
        
        std_threshold = label_config["std_threshold"]
        
        if "type_map" in selection_config:
            type_map = selection_config["type_map"]
        else:
            type_map = []
        
        op_out = OPIO(
            {
                "numb_iters": numb_iters,
                "numb_walkers": numb_walkers,
                "numb_models": numb_models,
                "confs": conf_list,
                "walker_tags": walker_tags,
                "model_tags": model_tags,

                "exploration_config": exploration_config,
                "cv_config": cv_config,
                "trust_lvl_1": trust_lvl_1_list,
                "trust_lvl_2": trust_lvl_2_list,
                "cluster_threshold": cluster_threshold_list,
                "angular_mask": angular_mask,
                "weights": weights,
                "numb_cluster_upper": selection_config.pop("numb_cluster_upper"),
                "numb_cluster_lower": selection_config.pop("numb_cluster_lower"),
                "max_selection": selection_config.pop("max_selection"),
                "numb_cluster_threshold": selection_config.pop("numb_cluster_threshold"),
                "std_threshold": std_threshold,
                "dt": dt,
                "output_freq": output_freq,
                "slice_mode": selection_config.pop("slice_mode"),
                "type_map": type_map,
                "label_config": label_config,
                "train_config": train_config
            }
        )
        return op_out