from dflow.python import (
    OP,
    OPIO,
    OPIOSign,
    Artifact
)

import os, sys, json, shutil, logging
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
from copy import deepcopy
from rid.task.builder import RestrainedMDTaskBuilder
from rid.utils import load_txt, load_json
from rid.constants import model_tag_fmt, init_conf_name, walker_tag_fmt


logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def prep_confs(confs, numb_walkers):
    numb_confs = len(confs)
    conf_list = []
    if numb_confs < numb_walkers:
        logger.info("Number of confs is smaller than number of walkers. Copy replicas up to number of walkers.")
        for idx in range(numb_walkers):
            shutil.copyfile(confs[idx%numb_confs], init_conf_name.format(idx=idx))
    elif numb_confs > numb_walkers:
        logger.info("Number of confs is greater than number of walkers. Only use the fist `numb_walkers` confs.")
        for idx in range(numb_walkers):
            shutil.copyfile(confs[idx], init_conf_name.format(idx=idx))
    else:
        for idx in range(numb_walkers):
            shutil.copyfile(confs[idx], init_conf_name.format(idx=idx))
    for idx in range(numb_walkers):
        conf_list.append(Path(init_conf_name.format(idx=idx)))
    return conf_list


class PrepRiD(OP):

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
                "confs": Artifact(List[Path]),
                "walker_tags": List,
                "model_tags": List,
                
                "exploration_gmx_config": Dict,
                "cv_config": Dict,
                "trust_lvl_1": List[float],
                "trust_lvl_2": List[float],
                "cluster_threshold": float,
                "angular_mask": List,
                "weights": List,
                "numb_cluster_upper": int,
                "numb_cluster_lower": int,
                "max_selection": int,
                "numb_cluster_threshold": int,
                "dt": float,
                "slice_mode": str,
                "label_gmx_config": Dict,
                "kappas": List,
                "train_config": Dict
            }
        )

    @OP.exec_sign_check
    def execute(
        self,
        op_in: OPIO,
    ) -> OPIO:
        jdata = deepcopy(load_json(op_in["rid_config"]))
        numb_walkers = jdata.pop("numb_walkers")
        train_config = jdata.pop("Train")
        numb_models = train_config.pop("numb_models")
        numb_iters = jdata.pop("numb_iters")
        conf_list = prep_confs(op_in["confs"], numb_walkers)

        walker_tags = []
        model_tags = []
        for idx in range(numb_walkers):
            walker_tags.append(walker_tag_fmt.format(idx=idx))
        for idx in range(numb_models):
            model_tags.append(model_tag_fmt.format(idx=idx))
        
        exploration_gmx_config = jdata.pop("ExploreMDConfig")
        dt = exploration_gmx_config["dt"]
        cv_config = jdata.pop("CV")
        angular_mask = cv_config.pop("angular_mask")
        weights = cv_config.pop("weights")
        
        selection_config = jdata.pop("SelectorConfig")

        label_md_config = jdata.pop("LabelMDConfig")
        kappas = label_md_config.pop("kappas")
        
        trust_lvl_1 = jdata.pop("trust_lvl_1")
        trust_lvl_2 = jdata.pop("trust_lvl_2")
        trust_lvl_1_list = [trust_lvl_1 for _ in range(numb_walkers)]
        trust_lvl_2_list = [trust_lvl_2 for _ in range(numb_walkers)]
        
        op_out = OPIO(
            {
                "numb_iters": numb_iters,
                "numb_walkers": numb_walkers,
                "numb_models": numb_models,
                "confs": conf_list,
                "walker_tags": walker_tags,
                "model_tags": model_tags,

                "exploration_gmx_config": exploration_gmx_config,
                "cv_config": cv_config,
                "trust_lvl_1": trust_lvl_1_list,
                "trust_lvl_2": trust_lvl_2_list,
                "cluster_threshold": selection_config.pop("cluster_threshold"),
                "angular_mask": angular_mask,
                "weights": weights,
                "numb_cluster_upper": selection_config.pop("numb_cluster_upper"),
                "numb_cluster_lower": selection_config.pop("numb_cluster_lower"),
                "max_selection": selection_config.pop("max_selection"),
                "numb_cluster_threshold": selection_config.pop("numb_cluster_threshold"),
                "dt": dt,
                "slice_mode": selection_config.pop("slice_mode"),
                "label_gmx_config": label_md_config,
                "kappas": kappas,
                "train_config": jdata.pop("Train")
            }
        )
        return op_out