import json
from pathlib import Path
from typing import List, Union, Optional

from dflow import (
    Workflow,
    Step,
    upload_artifact
)

from dflow.python import upload_packages
from rid import SRC_ROOT
upload_packages.append(SRC_ROOT)

from rid.utils import normalize_resources
from .submit import prep_rid_op


def resubmit_rid(
        workflow_id: str,
        confs: Union[str, List[str]],
        topology: Optional[str],
        rid_config: str,
        machine_config: str,
        models: Optional[Union[str, List[str]]] = None,
        index_file: Optional[str] = None,
        dp_files: Optional[List[str]] = None,
        inputfile: Optional[List[str]] = None,
        forcefield: Optional[str] = None,
        cv_file: Optional[List[str]] = None
    ):
    with open(machine_config, "r") as mcg:
        machine_config_dict = json.load(mcg)
    resources = machine_config_dict["resources"]
    tasks = machine_config_dict["tasks"]
    normalized_resources = {}
    for resource_type in resources.keys():
        normalized_resources[resource_type] = normalize_resources(resources[resource_type])

    rid_op = prep_rid_op(
        prep_exploration_config = normalized_resources[tasks["prep_exploration_config"]],
        run_exploration_config = normalized_resources[tasks["run_exploration_config"]],
        prep_label_config = normalized_resources[tasks["prep_label_config"]],
        run_label_config = normalized_resources[tasks["run_label_config"]],
        prep_select_config = normalized_resources[tasks["prep_select_config"]],
        run_select_config = normalized_resources[tasks["run_select_config"]],
        prep_data_config = normalized_resources[tasks["prep_data_config"]],
        run_train_config = normalized_resources[tasks["run_train_config"]],
        workflow_steps_config = normalized_resources[tasks["workflow_steps_config"]]
    )

    if isinstance(confs, str):
        confs_artifact = upload_artifact(Path(confs), archive=None)
    elif isinstance(confs, List):
        confs_artifact = upload_artifact([Path(p) for p in confs], archive=None)
    else:
        raise RuntimeError("Invalid type of `confs`.")
    
    if models is None:
        models_artifact = None
    elif isinstance(models, str):
        models_artifact = upload_artifact(Path(models), archive=None)
    elif isinstance(models, List):
        models_artifact = upload_artifact([Path(p) for p in models], archive=None)
    else:
        raise RuntimeError("Invalid type of `confs`.")
    
    if index_file is None:
        index_file_artifact = None
    else:
        index_file_artifact = upload_artifact(Path(index_file), archive=None)
        
    if inputfile is None:
        inputfile_artifact = None
    else:
        inputfile_artifact = upload_artifact([Path(p) for p in inputfile], archive=None)
        
    if dp_files is None:
        dp_files_artifact = None
    elif isinstance(dp_files, str):
        dp_files_artifact = upload_artifact(Path(dp_files), archive=None)
    elif isinstance(dp_files, List):
        dp_files_artifact = upload_artifact([Path(p) for p in dp_files], archive=None)
    else:
        raise RuntimeError("Invalid type of `dp_files`.")
    
    if cv_file is None:
        cv_file_artifact = None
    elif isinstance(cv_file, str):
        cv_file_artifact = upload_artifact(Path(cv_file), archive=None)
    elif isinstance(cv_file, List):
        cv_file_artifact = upload_artifact([Path(p) for p in cv_file], archive=None)
    else:
        raise RuntimeError("Invalid type of `cv_file`.")
    
    if forcefield is None:
        forcefield_artifact = None
    else:
        forcefield_artifact = upload_artifact(Path(forcefield), archive=None)
    
    if topology is None:
        top_artifact = None
    else:
        top_artifact = upload_artifact(Path(topology), archive=None)
    rid_config = upload_artifact(Path(rid_config), archive=None)

    rid_steps = Step(
        "rid-procedure",
        rid_op,
        artifacts={
            "topology": top_artifact,
            "confs": confs_artifact,
            "rid_config": rid_config,
            "models": models_artifact,
            "forcefield": forcefield_artifact,
            "index_file": index_file_artifact,
            "inputfile": inputfile_artifact,
            "dp_files": dp_files_artifact,
            "cv_file": cv_file_artifact
        },
        parameters={}
        )
    old_workflow = Workflow(id=workflow_id)
    all_steps = old_workflow.query_step()

    succeeded_steps = []
    for step in all_steps:
        if step["type"] == "Pod":
            if step["phase"] == "Succeeded":
                if step["key"] != "prepare-rid":
                    succeeded_steps.append(step)
    wf = Workflow("reinforced-dynamics", pod_gc_strategy="OnPodSuccess", parallelism=30)
    wf.add(rid_steps)
    wf.submit(reuse_step=succeeded_steps)
