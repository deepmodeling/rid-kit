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
    argo_sequence
)
from dflow.python import(
    PythonOPTemplate,
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Slices
)
from rid.utils import init_executor


class Exploration(Steps):
    def __init__(
        self,
        name: str,
        prep_op: OP,
        run_op: OP,
        prep_config: Dict,
        run_config: Dict,
        upload_python_package = None
    ):
        self._input_parameters = {
            "trust_lvl_1" : InputParameter(type=float, value=2.0),
            "trust_lvl_2": InputParameter(type=float, value=3.0),
            "gmx_config" : InputParameter(type=Dict),
            "cv_config" : InputParameter(type=Dict),
            "task_names" : InputParameter(type=str)
        }        
        self._input_artifacts = {
            "models" : InputArtifact(optional=True),
            "topology" : InputArtifact(),
            "confs" : InputArtifact(),
        }
        self._output_parameters = {
            "cv_dim": OutputParameter(type=int)
        }
        self._output_artifacts = {
            "plm_out": OutputArtifact(),
            "md_log": OutputArtifact(),
            "trajectory": OutputArtifact()
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

        self = _exploration(
            self, 
            prep_op,
            run_op,
            prep_config = prep_config,
            run_config = run_config,
            upload_python_package = upload_python_package,
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
        prep_exploration_op : OP,
        run_exploration_op : OP,
        prep_config : Dict,
        run_config : Dict,
        upload_python_package : str = None,
    ):
    prep_config = deepcopy(prep_config)
    run_config = deepcopy(run_config)
    prep_template_config = prep_config.pop('template_config')
    run_template_config = run_config.pop('template_config')
    prep_executor = init_executor(prep_config.pop('executor'))
    run_executor = init_executor(run_config.pop('executor'))
    prep_exploration = Step(
        'prep-exploration',
        template=PythonOPTemplate(
            prep_exploration_op,
            python_packages = upload_python_package,
            slices=Slices(
                input_parameter=["task_name"],
                input_artifact=["conf"],
                output_artifact=["task_path"]),
            **prep_template_config,
        ),
        parameters={
            "trust_lvl_1" : exploration_steps.inputs.parameters['trust_lvl_1'],
            "trust_lvl_2": exploration_steps.inputs.parameters['trust_lvl_2'],
            "gmx_config" : exploration_steps.inputs.parameters['gmx_config'],
            "cv_config" : exploration_steps.inputs.parameters['cv_config'],
            "task_name": exploration_steps.inputs.parameters['task_names']
        },
        artifacts={
            "models" : exploration_steps.inputs.artifacts['models'],
            "topology" :exploration_steps.inputs.artifacts['topology'],
            "conf" : exploration_steps.inputs.artifacts['confs']
        },
        key = 'prep-exploration',
        with_param=argo_range(argo_len(exploration_steps.inputs.parameters['task_names'])),
        executor = prep_executor,
        **prep_config,
    )
    exploration_steps.add(prep_exploration)

    run_exploration = Step(
        'run-exploration',
        template=PythonOPTemplate(
            run_exploration_op,
            python_packages = upload_python_package,
             slices=Slices(
                input_artifact=["task_path"],
                output_artifact=["plm_out", "trajectory", "md_log"]),
            **run_template_config,
        ),
        parameters={
            "gmx_config" : exploration_steps.inputs.parameters["gmx_config"]
        },
        artifacts={
            "task_path" : prep_exploration.outputs.artifacts["task_path"],
            "models" : exploration_steps.inputs.artifacts['models']
        },
        key = "run-exploration",
        executor = run_executor,
        with_param=argo_range(argo_len(exploration_steps.inputs.parameters['task_names'])),
        **run_config,
    )
    exploration_steps.add(run_exploration)

    exploration_steps.outputs.parameters["cv_dim"].value_from_parameter = prep_exploration.outputs.parameters["cv_dim"]
    exploration_steps.outputs.artifacts["plm_out"]._from = run_exploration.outputs.artifacts["plm_out"]
    exploration_steps.outputs.artifacts["md_log"]._from = run_exploration.outputs.artifacts["md_log"]
    exploration_steps.outputs.artifacts["trajectory"]._from = run_exploration.outputs.artifacts["trajectory"]
    
    return exploration_steps