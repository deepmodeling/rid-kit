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
    