from distutils.command.config import dump_file
import json
from pathlib import Path
from typing import List, Union, Optional
from rid.utils import load_json
from copy import deepcopy
import os

from dflow import (
    Workflow,
    Step,
    upload_artifact
)

from dflow.python import(
    PythonOPTemplate,
    upload_packages,
    OP
)
from rid import SRC_ROOT
upload_packages.append(SRC_ROOT)

from rid.utils import normalize_resources,init_executor
from rid.op.prep_rid import PrepRiD
from rid.superop.exploration import Exploration
from rid.op.prep_exploration import PrepExplore
from rid.op.run_exploration import RunExplore

def explore_rid(
        confs: Union[str, List[str]],
        topology: Optional[str],
        rid_config: str,
        machine_config: str,
        models: Optional[Union[str, List[str]]] = None,
        forcefield: Optional[str] = None,
        index_file: Optional[str] = None,
        dp_files: Optional[List[str]] = [],
        otherfiles: Optional[List[str]] = None
    ):
    with open(machine_config, "r") as mcg:
        machine_config_dict = json.load(mcg)
    resources = machine_config_dict["resources"]
    tasks = machine_config_dict["tasks"]
    normalized_resources = {}
    for resource_type in resources.keys():
        normalized_resources[resource_type] = normalize_resources(resources[resource_type])
        
    if isinstance(confs, str):
        confs_artifact = upload_artifact(Path(confs), archive=None)
    elif isinstance(confs, List):
        confs_artifact = upload_artifact([Path(p) for p in confs], archive=None)
    else:
        raise RuntimeError("Invalid type of `confs`.")

    if index_file is None:
        index_file_artifact = None
    else:
        index_file_artifact = upload_artifact(Path(index_file), archive=None)
    
    jdata = deepcopy(load_json(rid_config))
    cv_config = jdata["CV"]
    explore_config = jdata["ExploreMDConfig"]
    
    inputfiles = []
    if "inputfile" in explore_config:
        inputfiles.append(explore_config["inputfile"])

    fe_models = []
    assert isinstance(jdata["init_models"],list), "model input should be list."
    for model in jdata["init_models"]:
        fe_models.append(model)
    
    cvfiles = []
    if "cv_file" in jdata["CV"]:
        assert isinstance(jdata["CV"]["cv_file"],list), "CV file input should be list."
        for file in jdata["CV"]["cv_file"]:
            cvfiles.append(file)
            
    dp_models = []
    if "dp_model" in explore_config:
        assert isinstance(explore_config["dp_model"],list), "model input should be list."
        for model in explore_config["dp_model"]:
            dp_models.append(model)
                
    inputfile_list = []
    cvfile_list = []
    model_list = []
    dpfile_list = []
    if otherfiles is not None:
        for file in otherfiles:
            if os.path.basename(file) in inputfiles:
                inputfile_list.append(file)
            elif os.path.basename(file) in cvfiles:
                cvfile_list.append(file)
                
    if dp_files is not None:
        for dp_file in dp_files:
            dpfile_list.append(dp_file)
            
    if models is not None:
        for model in models:
            if os.path.basename(model) in fe_models:
                model_list.append(model)
            elif os.path.basename(model) in dp_models:
                dpfile_list.append(model)
            
    if len(inputfile_list) == 0:
        inputfile_artifact = None
    else:
        inputfile_artifact = upload_artifact([Path(p) for p in inputfile_list], archive=None)
        
    if len(model_list) == 0:
        models_artifact = None
    else:
        models_artifact = upload_artifact([Path(p) for p in model_list], archive=None)
        
    if len(cvfile_list) == 0:
        cv_file_artifact = None
    else:
        cv_file_artifact = upload_artifact([Path(p) for p in cvfile_list], archive=None)
        
    if len(dpfile_list) == 0:
        dp_files_artifact = None
    elif isinstance(dp_files, List):
        dp_files_artifact = upload_artifact([Path(p) for p in dpfile_list], archive=None)
    else:
        raise RuntimeError("Invalid type of `dp_files`.")
    
    if forcefield is None:
        forcefield_artifact = None
    else:
        forcefield_artifact = upload_artifact(Path(forcefield), archive=None)
        
    if topology is None:
        top_artifact = None
    else:
        top_artifact = upload_artifact(Path(topology), archive=None)
        
    rid_config = upload_artifact(Path(rid_config), archive=None)
        
    _step_config = normalized_resources[tasks["workflow_steps_config"]]
    _step_config = deepcopy(_step_config)
    step_template_config = _step_config.pop('template_config')
    step_executor = init_executor(_step_config.pop('executor'))
    
    prep_rid = Step(
        name = 'prepare-rid',
        template=PythonOPTemplate(
            PrepRiD,
            python_packages = None,
            **step_template_config,
        ),
        parameters={},
        artifacts={
            "confs":confs_artifact,
            "rid_config" : rid_config
        },
        key = 'prepare-rid',
        executor = step_executor,
        **_step_config,
    )

    explore_op = Exploration(
        "explore",
        PrepExplore,
        RunExplore,
        prep_config = normalized_resources[tasks["prep_exploration_config"]],
        run_config = normalized_resources[tasks["run_exploration_config"]],
        retry_times=1)


    exploration_steps = Step("rid-explore",
            explore_op,
            artifacts={
                "models": models_artifact,
                "forcefield": forcefield_artifact,
                "topology": top_artifact,
                "inputfile": inputfile_artifact,
                "confs": prep_rid.outputs.artifacts["confs"],
                "index_file": index_file_artifact,
                "dp_files": dp_files_artifact,
                "cv_file": cv_file_artifact
            },
            parameters={
                "trust_lvl_1" : prep_rid.outputs.parameters["trust_lvl_1"],
                "trust_lvl_2": prep_rid.outputs.parameters["trust_lvl_2"],
                "exploration_config" :explore_config,
                "cv_config": cv_config,
                "task_names": prep_rid.outputs.parameters["walker_tags"],
                "block_tag" : "iter-001"
            },
        )
    wf = Workflow("rid-exploration", pod_gc_strategy="OnPodSuccess", parallelism=50)
    wf.add(prep_rid)
    wf.add(exploration_steps)
    wf.submit()
