from typing import Dict, List
from copy import deepcopy
from dflow import (
    InputParameter,
    OutputParameter,
    Inputs,
    InputArtifact,
    Outputs,
    OutputArtifact,
    Step,
    Steps,
    argo_range,
    argo_len,
)
from dflow.python import(
    PythonOPTemplate,
    OP,
    Slices
)
from rid.utils import init_executor


class MCMC(Steps):
    
    r"""" MCMC SuperOP.
    This SuperOP combines MCMC_Run OP and MCMC_Plot OP.
    """
    def __init__(
        self,
        name: str,
        mcmc_run_op: OP,
        mcmc_plot_op: OP,
        run_config: Dict,
        plot_config: Dict,
        upload_python_package = None,
        retry_times = None
    ):
        self._input_parameters = {
            "mcmc_config" : InputParameter(type=Dict),
            "task_names" : InputParameter(type=List[str]),
            "block_tag" : InputParameter(type=str, value="")
        }        
        self._input_artifacts = {
            "models" : InputArtifact(),
            "plm_out": InputArtifact()
        }
        self._output_parameters = {}
        self._output_artifacts = {
            "mcmc_fig": OutputArtifact()
        }

        super().__init__(        
                name=name,
                inputs=Inputs(
                    parameters=self._input_parameters,
                    artifacts=self._input_artifacts
                ),
                outputs=Outputs(
                    parameters=self._output_parameters,
                    artifacts=self._output_artifacts
                ),
            )

        step_keys = {
            "mcmc_run": "{}-mcmc-run".format(self.inputs.parameters["block_tag"]),
            "mcmc_plot": "{}-mcmc-plot".format(self.inputs.parameters["block_tag"]),
        }

        self = _mcmc(
            self, 
            step_keys,
            mcmc_run_op,
            mcmc_plot_op,
            run_config = run_config,
            plot_config = plot_config,
            upload_python_package = upload_python_package,
            retry_times = retry_times
        )            
    
    @property
    def input_parameters(self):
        return self._input_parameters

    @property
    def input_artifacts(self):
        return self._input_artifacts

    @property
    def output_parameters(self):
        return self._output_parameters

    @property
    def output_artifacts(self):
        return self._output_artifacts

    @property
    def keys(self):
        return self._keys


def _mcmc(
        mcmc_steps,
        step_keys,
        mcmc_run_op : OP,
        mcmc_plot_op : OP,
        run_config : Dict,
        plot_config : Dict,
        upload_python_package : str = None,
        retry_times: int = None
    ):
    run_config = deepcopy(run_config)
    plot_config = deepcopy(plot_config)
    
    run_template_config = run_config.pop('template_config')
    plot_template_config = plot_config.pop('template_config')
    run_executor = init_executor(run_config.pop('executor'))
    plot_executor = init_executor(plot_config.pop('executor'))
    
    run_merge = False
    if run_executor is not None:
        run_merge = run_executor.merge_sliced_step
    if run_merge:
        mcmc_run = Step(
        'mcmc-run',
        template=PythonOPTemplate(
            mcmc_run_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            slices=Slices("{{item}}",
                input_parameter=["task_names"],
                input_artifact=["models"],
                output_artifact=["mcmc_1cv", "mcmc_2cv"]
            ),
            **run_template_config,
        ),
        parameters={
            "mcmc_config" : mcmc_steps.inputs.parameters['mcmc_config'],
            "task_names": mcmc_steps.inputs.parameters['task_names']
        },
        artifacts={
            "models" : mcmc_steps.inputs.artifacts['models']
        },
        key = step_keys["mcmc_run"]+"-{{item}}",
        with_param=argo_range(argo_len(mcmc_steps.inputs.parameters['task_names'])),
        executor = run_executor,
        **run_config
    )
    else:
        mcmc_run = Step(
            'mcmc-run',
            template=PythonOPTemplate(
                mcmc_run_op,
                python_packages = upload_python_package,
                retry_on_transient_error = retry_times,
                slices=Slices(sub_path = True,
                    input_parameter=["task_names"],
                    input_artifact=["models"],
                    output_artifact=["mcmc_1cv", "mcmc_2cv"]
                ),
                **run_template_config,
            ),
            parameters={
                "mcmc_config" : mcmc_steps.inputs.parameters['mcmc_config'],
                "task_names": mcmc_steps.inputs.parameters['task_names']
            },
            artifacts={
                "models" : mcmc_steps.inputs.artifacts['models']
            },
            key = step_keys["mcmc_run"]+"-{{item.order}}",
            executor = run_executor,
            **run_config
    )
    mcmc_steps.add(mcmc_run)

    mcmc_plot = Step(
        'mcmc-plot',
        template=PythonOPTemplate(
            mcmc_plot_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            **plot_template_config,
        ),
        parameters={
            "mcmc_config" : mcmc_steps.inputs.parameters["mcmc_config"]
        },
        artifacts={
            "mcmc_1cv": mcmc_run.outputs.artifacts['mcmc_1cv'],
            "mcmc_2cv": mcmc_run.outputs.artifacts['mcmc_2cv'],
            "plm_out":  mcmc_steps.inputs.artifacts['plm_out']
        },
        key = step_keys["mcmc_plot"],
        executor = plot_executor,
        **run_config
    )

    mcmc_steps.add(mcmc_plot)

    mcmc_steps.outputs.artifacts["mcmc_fig"]._from = mcmc_plot.outputs.artifacts["mcmc_fig"]
    
    return mcmc_steps