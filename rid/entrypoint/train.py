import json
from pathlib import Path
from rid.utils import load_json
from copy import deepcopy

from dflow import (
    Workflow,
    Step,
    Steps,
    upload_artifact,
    argo_range,
    argo_len
)

from dflow.python import(
    PythonOPTemplate,
    upload_packages,
    Slices,
)

from rid import SRC_ROOT
upload_packages.append(SRC_ROOT)

from rid.utils import normalize_resources,init_executor
from rid.op.run_train import TrainModel
from rid.constants import model_tag_fmt


def train_rid(
        data: str,
        rid_config: str,
        machine_config: str
    ):
    with open(machine_config, "r") as mcg:
        machine_config_dict = json.load(mcg)
    resources = machine_config_dict["resources"]
    tasks = machine_config_dict["tasks"]
    normalized_resources = {}
    for resource_type in resources.keys():
        normalized_resources[resource_type] = normalize_resources(resources[resource_type])

    jdata = deepcopy(load_json(rid_config))
    cv_config = jdata["CV"]
    train_config = jdata["Train"]
    numb_models = train_config["numb_models"]
    angular_mask = cv_config["angular_mask"]
    
    run_train_config = normalized_resources[tasks["run_train_config"]]
    run_train_config = deepcopy(run_train_config)
    train_template_config = run_train_config.pop('template_config')
    train_executor = init_executor(run_train_config.pop('executor'))
    
    model_tags = []
    for idx in range(numb_models):
        model_tags.append(model_tag_fmt.format(idx=idx))
        
    data_artifact = upload_artifact(Path(data), archive=None)
    
    train_op = TrainModel
    train_steps=Step(
        "train",
        template=PythonOPTemplate(
            train_op,
            python_packages = None,
            retry_on_transient_error = 1,
            slices=Slices("{{item}}",
                input_parameter=["model_tag"],
                output_artifact=["model","train_fig"]),
            **train_template_config,
        ),
        parameters={
            "model_tag": model_tags,
            "angular_mask": angular_mask,
            "train_config": train_config,
        },
        artifacts={
            "data": data_artifact
        },
        executor = train_executor,
        with_param = argo_range(len(model_tags)),
        key = "run-train"+"-{{item}}",
        **run_train_config,
    )
    wf = Workflow("rid-train", pod_gc_strategy="OnPodSuccess", parallelism=10)
    wf.add(train_steps)
    wf.submit()
