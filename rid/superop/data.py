from typing import Dict
from copy import deepcopy
from dflow import (
    InputParameter,
    Inputs,
    InputArtifact,
    Outputs,
    OutputArtifact,
    Step,
    Steps
)
from dflow.python import(
    PythonOPTemplate,
    OP
)
from rid.utils import init_executor


class DataGenerator(Steps):
    
    r""" Date generator SuperOP.
    This SuperpOP combines CollectData OP and MergeData OP to process data for training.
    """
    def __init__(
        self,
        name: str,
        collect_op: OP,
        merge_op: OP,
        run_config: Dict,
        upload_python_package = None,
        retry_times = None
    ):
        self._input_parameters = {
            "block_tag" : InputParameter(type=str, value="")
        }        
        self._input_artifacts = {
            "cv_forces": InputArtifact(),
            "data_old": InputArtifact(optional=True)
        }
        self._output_parameters = {}
        self._output_artifacts = {
            "data": OutputArtifact()
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

        self = _gen_data(
            self, 
            collect_op,
            merge_op,
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


def _gen_data(
        data_steps,
        collect_op : OP,
        merge_op : OP,
        run_config : Dict,
        upload_python_package : str = None,
        retry_times: int = None
    ):
    run_config = deepcopy(run_config)
    run_template_config = run_config.pop('template_config')
    run_executor = init_executor(run_config.pop('executor'))

    collect_data = Step(
        'collect-data',
        template=PythonOPTemplate(
            collect_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            **run_template_config,
        ),
        parameters={},
        artifacts={
            "cv_forces": data_steps.inputs.artifacts['cv_forces'],
        },
        key = '{}-collect-data'.format(data_steps.inputs.parameters["block_tag"]),
        executor = run_executor,
        **run_config,
    )
    data_steps.add(collect_data)
    
    merge_data = Step(
        'merge-data',
        template=PythonOPTemplate(
            merge_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            **run_template_config,
        ),
        parameters={},
        artifacts={
            "data_old": data_steps.inputs.artifacts["data_old"],
            "data_new": collect_data.outputs.artifacts["data_new"]
        },
        key = '{}-merge-data'.format(data_steps.inputs.parameters["block_tag"]),
        executor = run_executor,
        **run_config,
    )
    data_steps.add(merge_data)

    data_steps.outputs.artifacts["data"]._from = merge_data.outputs.artifacts["data_raw"]
   
    return data_steps
