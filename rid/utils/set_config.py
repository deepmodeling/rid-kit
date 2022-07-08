from dflow.plugins.lebesgue import LebesgueExecutor


def init_executor(
        executor_dict,
):
    if executor_dict is None:
        return None
    etype = executor_dict.pop('type')
    if etype == "lebesgue_v2":
        return LebesgueExecutor(**executor_dict)
    else:
        raise RuntimeError('unknown executor type', etype)    