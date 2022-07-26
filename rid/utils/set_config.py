from typing import Dict, List
from dflow.plugins.lebesgue import LebesgueExecutor
from dflow import SlurmRemoteExecutor


def init_executor(
        executor_dict,
    ):
    if executor_dict is None:
        return None
    etype = executor_dict.pop('type')
    if etype.lower() == "lebesgue_v2":
        return LebesgueExecutor(**executor_dict)
    elif etype.lower() == "slurm":
        return SlurmRemoteExecutor(**executor_dict)
    else:
        raise RuntimeError('unknown executor type', etype)    


def normalize_resources(config_dict: Dict):
    template_dict = {}
    template_dict["template_config"] = config_dict.get("template_config", {})
    template_dict["executor"] = config_dict.get("executor", None)
    if template_dict["executor"] is None:
        assert ("image" in template_dict["template_config"].keys()) and \
            template_dict["template_config"]["image"] is not None
    if template_dict["executor"] is not None and template_dict["executor"]["type"] == "slurm":
        header_list = template_dict["executor"]["header"]
        template_dict["executor"]["header"] = "\n".join(header_list)
    return template_dict

