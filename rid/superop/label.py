from typing import Dict, List
from copy import deepcopy
from dflow import (
    InputParameter,
    OutputParameter,
    Inputs,
    InputArtifact,
    Outputs,
    OutputArtifact,
    Workflow,
    Step,
    Steps,
    upload_artifact,
    download_artifact,
    argo_range,
    argo_len,
    argo_sequence,
)
from dflow.python import(
    PythonOPTemplate,
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Slices,
)
from rid.utils import init_executor


class Label(Steps):
    def __init__(
        self,
        name: str,
        prep_op: OP,
        run_op: OP,
        post_op: OP,
        prep_config: Dict,
        run_config: Dict,
        upload_python_package = None
    ):

        self._input_parameters = {
            "gmx_config": InputParameter(type=Dict),
            "cv_config": InputParameter(type=Dict),
            "kappas": InputParameter(type=List[float]),
            "angular_mask": InputParameter(type=List),
            "tail": InputParameter(type=float, value=0.9),
            "task_names": InputParameter(type=str),
        }        
        self._input_artifacts = {
            "topology": InputArtifact(),
            "confs": InputArtifact(),
            "at": InputArtifact()
        }
        self._output_parameters = {
        }
        self._output_artifacts = {
            "md_log": OutputArtifact(),
            "forces": OutputArtifact()
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

        self = _label(
            self, 
            prep_op,
            run_op,
            post_op,
            prep_config = prep_config,
            run_config = run_config,
            post_config = prep_config,
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


def _label(
        label_steps,
        prep_label_op : OP,
        run_label_op : OP,
        post_label_op : OP,
        prep_config : Dict,
        run_config : Dict,
        post_config : Dict,
        upload_python_package : str = None,
    ):
    prep_config = deepcopy(prep_config)
    run_config = deepcopy(run_config)
    post_config = deepcopy(post_config)
    prep_template_config = prep_config.pop('template_config')
    run_template_config = run_config.pop('template_config')
    post_template_config = post_config.pop('template_config')
    prep_executor = init_executor(prep_config.pop('executor'))
    run_executor = init_executor(run_config.pop('executor'))
    post_executor = init_executor(post_config.pop('executor'))

    prep_label = Step(
        'prep-label',
        template=PythonOPTemplate(
            prep_label_op,
            python_packages = upload_python_package,
            slices=Slices(
                input_parameter=["task_name"],
                input_artifact=["conf", "at"],
                output_artifact=["task_path"]),
            **prep_template_config,
        ),
        parameters={
            "gmx_config": label_steps.inputs.parameters['gmx_config'],
            "cv_config": label_steps.inputs.parameters['cv_config'],
            "task_name": label_steps.inputs.parameters['task_names'],
            "kappas": label_steps.inputs.parameters['kappas']
        },
        artifacts={
            "topology": label_steps.inputs.artifacts['topology'],
            "conf": label_steps.inputs.artifacts['confs'],
            "at": label_steps.inputs.artifacts['at']
        },
        key = 'prep-label',
        executor = prep_executor,
        with_param=argo_range(argo_len(label_steps.inputs.parameters['task_names'])),
        **prep_config,
    )
    label_steps.add(prep_label)

    run_label = Step(
        'run-label',
        template=PythonOPTemplate(
            run_label_op,
            python_packages = upload_python_package,
            slices=Slices(
                input_artifact=["task_path"],
                output_artifact=["plm_out", "md_log"]),
            **run_template_config,
        ),
        parameters={
            "gmx_config": label_steps.inputs.parameters["gmx_config"]
        },
        artifacts={
            "task_path": prep_label.outputs.artifacts["task_path"]
        },
        key = "run-label",
        executor = run_executor,
        with_param=argo_range(argo_len(label_steps.inputs.parameters['task_names'])),
        **run_config,
    )
    label_steps.add(run_label)

    post_label = Step(
        'cpmf',
        template=PythonOPTemplate(
            post_label_op,
            python_packages = upload_python_package,
            slices=Slices(
                input_parameter=["task_name"],
                input_artifact=["plm_out", "at"],
                output_artifact=["forces"]),
            **run_template_config,
        ),
        parameters={
            "task_name": label_steps.inputs.parameters['task_names'],
            "kappas": label_steps.inputs.parameters['kappas'],
            "tail": label_steps.inputs.parameters['tail'],
            "angular_mask": label_steps.inputs.parameters['angular_mask']
        },
        artifacts={
            "plm_out": run_label.outputs.artifacts["plm_out"],
            "at": label_steps.inputs.artifacts['at']
        },
        key = "cpmf",
        executor = run_executor,
        with_param=argo_range(argo_len(label_steps.inputs.parameters['task_names'])),
        **run_config,
    )
    label_steps.add(post_label)

    label_steps.outputs.artifacts["forces"]._from = post_label.outputs.artifacts["forces"]
    label_steps.outputs.artifacts["md_log"]._from = run_label.outputs.artifacts["md_log"]
    
    return label_steps
