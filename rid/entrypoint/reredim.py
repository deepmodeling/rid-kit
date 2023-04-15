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

from dflow.python import upload_packages
from rid import SRC_ROOT
upload_packages.append(SRC_ROOT)

from rid.utils import normalize_resources
from rid.superop.mcmc import MCMC
from rid.op.mcmc_run import MCMCRun
from rid.op.mcmc_plot import MCMCPlot

def reredim_rid(
        workflow_id: str,
        rid_config: str,
        machine_config: str,
        models: Optional[Union[str, List[str]]] = None,
        pod: Optional[str] = None
    ):
    with open(machine_config, "r") as mcg:
        machine_config_dict = json.load(mcg)
    resources = machine_config_dict["resources"]
    tasks = machine_config_dict["tasks"]
    normalized_resources = {}
    for resource_type in resources.keys():
        normalized_resources[resource_type] = normalize_resources(resources[resource_type])

    mcmc_op = MCMC(
        "mcmc",
        MCMCRun,
        MCMCPlot,
        run_config = normalized_resources[tasks["mcmc_run_config"]],
        plot_config = normalized_resources[tasks["mcmc_plot_config"]],
        retry_times=None)
    
    jdata = deepcopy(load_json(rid_config))
    mcmc_config = jdata["MCMC_Config"]
    
               
    fe_models = []
    assert isinstance(jdata["init_models"],list), "model input should be list."
    for model in jdata["init_models"]:
        fe_models.append(model)
                
    model_list = []
            
    if models is not None:
        for model in models:
            if os.path.basename(model) in fe_models:
                model_list.append(model)

    if len(model_list) == 0:
        models_artifact = None
    else:
        models_artifact = upload_artifact([Path(p) for p in model_list], archive=None)
    
    task_names = []
    for index in range(len(model_list)):
        task_names.append("%03d"%index)

    rid_steps = Step("rid-mcmc",
            mcmc_op,
            artifacts={
                "models": models_artifact
            },
            parameters={
                "mcmc_config": mcmc_config,
                "task_names": task_names,
                "block_tag" : "000"
            },
        )
    
    old_workflow = Workflow(id=workflow_id)
    all_steps = old_workflow.query_step()

    succeeded_steps = []
    restart_flag = 1
    for step in all_steps:
        if step["type"] == "Pod":
            pod_key = step["key"]
            if pod_key is not None:
                pod_key_list = pod_key.split("-")
                pod_step = "-".join(pod_key_list[1:-1])
                if pod is not None:
                    if pod_step == pod:
                        restart_flag = 0
                else:
                    if step["phase"] != "Succeeded":
                        restart_flag = 0
                    else:
                        restart_flag = 1
            
            if restart_flag == 1:
                succeeded_steps.append(step)
                
    wf = Workflow("rid-mcmc-continue", pod_gc_strategy="OnPodSuccess", parallelism=10)
    wf.add(rid_steps)
    wf.submit(reuse_step=succeeded_steps)
