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


class Exploration(Steps):
    
    r"""" Exploration SuperOP.
    This SuperOP combines PrepExplore OP and RunExplore OP.
    """
    def __init__(
        self,
        name: str,
        prep_op: OP,
        run_op: OP,
        prep_config: Dict,
        run_config: Dict,
        upload_python_package = None,
        retry_times = None
    ):
        self._input_parameters = {
            "trust_lvl_1" : InputParameter(type=List[float], value=2.0),
            "trust_lvl_2": InputParameter(type=List[float], value=3.0),
            "exploration_config" : InputParameter(type=Dict),
            "cv_config" : InputParameter(type=Dict),
            "task_names" : InputParameter(type=List[str]),
            "block_tag" : InputParameter(type=str, value="")
        }        
        self._input_artifacts = {
            "models" : InputArtifact(optional=True),
            "forcefield": InputArtifact(optional=True),
            "topology" : InputArtifact(optional=True),
            "inputfile": InputArtifact(optional=True),
            "confs" : InputArtifact(),
            "index_file": InputArtifact(optional=True),
            "dp_files": InputArtifact(optional=True),
            "cv_file": InputArtifact(optional=True)
        }
        self._output_parameters = {
            "cv_dim": OutputParameter(type=List[int])
        }
        self._output_artifacts = {
            "plm_out": OutputArtifact(),
            "md_log": OutputArtifact(),
            "trajectory": OutputArtifact(),
            "conf_outs": OutputArtifact()
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
            "prep_exploration": "{}-prep-exploration".format(self.inputs.parameters["block_tag"]),
            "run_exploration": "{}-run-exploration".format(self.inputs.parameters["block_tag"]),
        }

        self = _exploration(
            self, 
            step_keys,
            prep_op,
            run_op,
            prep_config = prep_config,
            run_config = run_config,
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


def _exploration(
        exploration_steps,
        step_keys,
        prep_exploration_op : OP,
        run_exploration_op : OP,
        prep_config : Dict,
        run_config : Dict,
        upload_python_package : str = None,
        retry_times: int = None
    ):
    prep_config = deepcopy(prep_config)
    run_config = deepcopy(run_config)
    prep_template_config = prep_config.pop('template_config')
    run_template_config = run_config.pop('template_config')
    prep_executor = init_executor(prep_config.pop('executor'))
    run_executor = init_executor(run_config.pop('executor'))
    prep_merge = False
    if prep_executor is not None:
        prep_merge = prep_executor.merge_sliced_step
    if prep_merge:
        prep_exploration = Step(
        'prep-exploration',
        template=PythonOPTemplate(
            prep_exploration_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            slices=Slices("{{item}}",
                input_parameter=["task_name", "trust_lvl_1", "trust_lvl_2"],
                input_artifact=["conf"],
                output_artifact=["task_path"],
                output_parameter=["cv_dim"]
            ),
            **prep_template_config,
        ),
        parameters={
            "trust_lvl_1" : exploration_steps.inputs.parameters['trust_lvl_1'],
            "trust_lvl_2": exploration_steps.inputs.parameters['trust_lvl_2'],
            "exploration_config" : exploration_steps.inputs.parameters['exploration_config'],
            "cv_config" : exploration_steps.inputs.parameters['cv_config'],
            "task_name": exploration_steps.inputs.parameters['task_names'],
            "block_tag": exploration_steps.inputs.parameters["block_tag"]
        },
        artifacts={
            "models" : exploration_steps.inputs.artifacts['models'],
            "topology" :exploration_steps.inputs.artifacts['topology'],
            "conf" : exploration_steps.inputs.artifacts['confs'],
            "cv_file": exploration_steps.inputs.artifacts['cv_file']
        },
        key = step_keys["prep_exploration"]+"-{{item}}",
        with_param=argo_range(argo_len(exploration_steps.inputs.parameters['task_names'])),
        executor = prep_executor,
        **prep_config
    )
    else:
        prep_exploration = Step(
            'prep-exploration',
            template=PythonOPTemplate(
                prep_exploration_op,
                python_packages = upload_python_package,
                retry_on_transient_error = retry_times,
                slices=Slices(sub_path = True,
                    input_parameter=["task_name", "trust_lvl_1", "trust_lvl_2"],
                    input_artifact=["conf"],
                    output_artifact=["task_path"],
                    output_parameter=["cv_dim"]
                ),
                **prep_template_config,
            ),
            parameters={
                "trust_lvl_1" : exploration_steps.inputs.parameters['trust_lvl_1'],
                "trust_lvl_2": exploration_steps.inputs.parameters['trust_lvl_2'],
                "exploration_config" : exploration_steps.inputs.parameters['exploration_config'],
                "cv_config" : exploration_steps.inputs.parameters['cv_config'],
                "task_name": exploration_steps.inputs.parameters['task_names'],
                "block_tag": exploration_steps.inputs.parameters["block_tag"]
            },
            artifacts={
                "models" : exploration_steps.inputs.artifacts['models'],
                "topology" :exploration_steps.inputs.artifacts['topology'],
                "conf" : exploration_steps.inputs.artifacts['confs'],
                "cv_file": exploration_steps.inputs.artifacts['cv_file']
            },
            key = step_keys["prep_exploration"]+"-{{item.order}}",
            executor = prep_executor,
            **prep_config
    )
    exploration_steps.add(prep_exploration)

    run_merge = False
    if run_executor is not None:
        run_merge = run_executor.merge_sliced_step
    if run_merge:
        run_exploration = Step(
        'run-exploration',
        template=PythonOPTemplate(
            run_exploration_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            slices=Slices("{{item}}",
                input_artifact=["task_path"],
                output_artifact=["plm_out", "bias_fig","model_devi_fig", "dp_model_devi_fig", "dp_model_devi", "dp_selected_indices","dp_selected_confs","projected_fig","trajectory", "md_log", "conf_out"]
            ),
            **run_template_config,
        ),
        parameters={
            "exploration_config" : exploration_steps.inputs.parameters["exploration_config"]
        },
        artifacts={
            "task_path" : prep_exploration.outputs.artifacts["task_path"],
            "forcefield": exploration_steps.inputs.artifacts['forcefield'],
            "models" : exploration_steps.inputs.artifacts['models'],
            "index_file": exploration_steps.inputs.artifacts['index_file'],
            "dp_files": exploration_steps.inputs.artifacts['dp_files'],
            "cv_file": exploration_steps.inputs.artifacts['cv_file'],
            "inputfile": exploration_steps.inputs.artifacts['inputfile']
        },
        key = step_keys["run_exploration"]+"-{{item}}",
        executor = run_executor,
        with_param=argo_range(argo_len(exploration_steps.inputs.parameters['task_names'])),
        **run_config
    )
    else:
        run_exploration = Step(
            'run-exploration',
            template=PythonOPTemplate(
                run_exploration_op,
                python_packages = upload_python_package,
                retry_on_transient_error = retry_times,
                slices=Slices(sub_path = True,
                    input_artifact=["task_path"],
                    output_artifact=["plm_out", "bias_fig","model_devi_fig","dp_model_devi_fig", "dp_model_devi", "dp_selected_indices","dp_selected_confs","projected_fig","trajectory", "md_log", "conf_out"]
                ),
                **run_template_config,
            ),
            parameters={
                "exploration_config" : exploration_steps.inputs.parameters["exploration_config"]
            },
            artifacts={
                "task_path" : prep_exploration.outputs.artifacts["task_path"],
                "forcefield": exploration_steps.inputs.artifacts['forcefield'],
                "models" : exploration_steps.inputs.artifacts['models'],
                "index_file": exploration_steps.inputs.artifacts['index_file'],
                "dp_files": exploration_steps.inputs.artifacts['dp_files'],
                "cv_file": exploration_steps.inputs.artifacts['cv_file'],
                "inputfile": exploration_steps.inputs.artifacts['inputfile']
            },
            key = step_keys["run_exploration"]+"-{{item.order}}",
            executor = run_executor,
            **run_config,
        )
    exploration_steps.add(run_exploration)

    exploration_steps.outputs.parameters["cv_dim"].value_from_parameter = prep_exploration.outputs.parameters["cv_dim"]
    exploration_steps.outputs.artifacts["plm_out"]._from = run_exploration.outputs.artifacts["plm_out"]
    exploration_steps.outputs.artifacts["md_log"]._from = run_exploration.outputs.artifacts["md_log"]
    exploration_steps.outputs.artifacts["trajectory"]._from = run_exploration.outputs.artifacts["trajectory"]
    exploration_steps.outputs.artifacts["conf_outs"]._from = run_exploration.outputs.artifacts["conf_out"]
    
    return exploration_steps