from typing import Dict
from dflow.plugins.lebesgue import LebesgueExecutor
from dflow.plugins.dispatcher import DispatcherExecutor
from dflow import SlurmRemoteExecutor
import os


def init_executor(
        executor_dict,
    ):
    if executor_dict is None:
        return None
    if "type" in executor_dict:
        etype = executor_dict.pop('type')
        if etype.lower() == "lebesgue_v2":
            return LebesgueExecutor(**executor_dict)
        elif etype.lower() == "slurm":
            return SlurmRemoteExecutor(**executor_dict)
        else:
            raise RuntimeError('unknown executor type', etype)    
    elif "machine_dict" in executor_dict:
        return DispatcherExecutor(**executor_dict)
    else:
        raise RuntimeError('unknown executor dict')   


def normalize_resources(config_dict: Dict):
    template_dict = {}
    template_dict["template_config"] = config_dict.get("template_config", {})
    template_dict["executor"] = config_dict.get("executor", None)
    if template_dict["executor"] is None:
        assert ("image" in template_dict["template_config"].keys()) and \
            template_dict["template_config"]["image"] is not None
    elif template_dict["executor"] is not None and not "type" in template_dict["executor"]:
        return template_dict
    elif template_dict["executor"] is not None and template_dict["executor"]["type"] == "slurm":
        header_list = template_dict["executor"]["header"]
        template_dict["executor"]["header"] = "\n".join(header_list)
    return template_dict