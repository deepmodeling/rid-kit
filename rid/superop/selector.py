from typing import Dict, List, Optional, Union
from copy import deepcopy
import numpy as np
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


class Selector(Steps):
    
    r""" Selector SuperOP.
    This SuperOP combines PrepSelect OP and RunSelect OP.    
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
            "label_config": InputParameter(type=Dict),
            "trust_lvl_1" : InputParameter(type=List[float], value=2.0),
            "trust_lvl_2": InputParameter(type=List[float], value=3.0),
            "cluster_threshold": InputParameter(type=List[float], value=1.0),
            "angular_mask": InputParameter(type=Optional[Union[np.ndarray, List]]),
            "weights": InputParameter(type=Optional[Union[np.ndarray, List]]),
            "numb_cluster_upper": InputParameter(type=Optional[float], value=None),
            "numb_cluster_lower": InputParameter(type=Optional[float], value=None),
            "max_selection": InputParameter(type=int),
            "dt": InputParameter(type=float, value=0.02),
            "output_freq": InputParameter(type=float, value=2500),
            "slice_mode": InputParameter(type=str, value="gmx"),
            "type_map": InputParameter(type=List, value=[]),
            "if_make_threshold": InputParameter(type=bool, value=False),
            "task_names" : InputParameter(type=List[str]),
            "block_tag" : InputParameter(type=str, value="")
        }        
        self._input_artifacts = {
            "models" : InputArtifact(optional=True),
            "plm_out": InputArtifact(),
            "xtc_traj": InputArtifact(),
            "topology": InputArtifact()
        }
        self._output_parameters = {
            "cluster_threshold": OutputParameter(type=List[int]),
            "numb_cluster": OutputParameter(type=List[int])
        }
        self._output_artifacts = {
            "cluster_selection_index": OutputArtifact(),
            "selected_confs": OutputArtifact(),
            "selected_cv_init": OutputArtifact(),
            "model_devi": OutputArtifact(),
            "selected_indices": OutputArtifact(),
            "selected_conf_tags": OutputArtifact()
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
            "prep_select": "{}-prep-select".format(self.inputs.parameters["block_tag"]),
            "run_select": "{}-run-select".format(self.inputs.parameters["block_tag"]),
            "post_select": "{}-post-select".format(self.inputs.parameters["block_tag"])
        }

        self = _select(
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


def _select(
        select_steps,
        step_keys,
        prep_select_op : OP,
        run_select_op : OP,
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
        prep_select = Step(
        'prep-select',
        template=PythonOPTemplate(
            prep_select_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            slices=Slices(
                "{{item}}",
                input_parameter=["cluster_threshold", "task_name"],
                input_artifact=["plm_out"],
                output_artifact=["cluster_fig","cluster_selection_index", "cluster_selection_data"],
                output_parameter=["cluster_threshold", "numb_cluster"]
                ),
                
            **prep_template_config,
        ),
        parameters={
            "cluster_threshold": select_steps.inputs.parameters['cluster_threshold'],
            "angular_mask": select_steps.inputs.parameters['angular_mask'],
            "weights": select_steps.inputs.parameters['weights'],
            "numb_cluster_upper": select_steps.inputs.parameters['numb_cluster_upper'],
            "numb_cluster_lower": select_steps.inputs.parameters['numb_cluster_lower'],
            "max_selection": select_steps.inputs.parameters['max_selection'],
            "if_make_threshold": select_steps.inputs.parameters['if_make_threshold'],
            "task_name": select_steps.inputs.parameters['task_names']
        },
        artifacts={
            "plm_out": select_steps.inputs.artifacts['plm_out']
        },
        key = step_keys["prep_select"]+"-{{item}}",
        executor = prep_executor,
        with_param=argo_range(argo_len(select_steps.inputs.parameters['task_names'])),
        **prep_config
    )  
    else:
        prep_select = Step(
        'prep-select',
        template=PythonOPTemplate(
            prep_select_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            slices=Slices(sub_path = True,
                input_parameter=["cluster_threshold", "task_name"],
                input_artifact=["plm_out"],
                output_artifact=["cluster_fig","cluster_selection_index", "cluster_selection_data"],
                output_parameter=["cluster_threshold", "numb_cluster"]
                ),
                
            **prep_template_config,
        ),
        parameters={
            "cluster_threshold": select_steps.inputs.parameters['cluster_threshold'],
            "angular_mask": select_steps.inputs.parameters['angular_mask'],
            "weights": select_steps.inputs.parameters['weights'],
            "numb_cluster_upper": select_steps.inputs.parameters['numb_cluster_upper'],
            "numb_cluster_lower": select_steps.inputs.parameters['numb_cluster_lower'],
            "max_selection": select_steps.inputs.parameters['max_selection'],
            "if_make_threshold": select_steps.inputs.parameters['if_make_threshold'],
            "task_name": select_steps.inputs.parameters['task_names']
        },
        artifacts={
            "plm_out": select_steps.inputs.artifacts['plm_out']
        },
        key = step_keys["prep_select"]+"-{{item.order}}",
        executor = prep_executor,
        **prep_config,
    )
    select_steps.add(prep_select)

    run_merge = False
    if run_executor is not None:
        run_merge = run_executor.merge_sliced_step
    if run_merge:
        run_select = Step(
        'run-select',
        template=PythonOPTemplate(
            run_select_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            slices=Slices(
                "int({{item}})",
                input_parameter=["task_name", "trust_lvl_1", "trust_lvl_2"],
                input_artifact=["cluster_selection_index", "cluster_selection_data", "xtc_traj", "topology"],
                output_artifact=["selected_confs", "selected_cv_init", "model_devi", "selected_indices","selected_conf_tags"]
            ),
            **run_template_config,
        ),
        parameters={
            "label_config":  select_steps.inputs.parameters["label_config"],
            "trust_lvl_1": select_steps.inputs.parameters["trust_lvl_1"],
            "trust_lvl_2": select_steps.inputs.parameters["trust_lvl_2"],
            "dt": select_steps.inputs.parameters["dt"],
            "output_freq": select_steps.inputs.parameters["output_freq"],
            "slice_mode": select_steps.inputs.parameters["slice_mode"],
            "type_map": select_steps.inputs.parameters["type_map"],
            "task_name": select_steps.inputs.parameters['task_names']
        },
        artifacts={
            "cluster_selection_index": prep_select.outputs.artifacts["cluster_selection_index"],
            "cluster_selection_data": prep_select.outputs.artifacts["cluster_selection_data"],
            "models": select_steps.inputs.artifacts["models"],
            "xtc_traj": select_steps.inputs.artifacts["xtc_traj"],
            "topology": select_steps.inputs.artifacts["topology"]
        },
        key = step_keys["run_select"]+"-{{item}}",
        executor = run_executor,
        with_param=argo_range(argo_len(select_steps.inputs.parameters["task_names"])),
        **run_config
    )
    else:
        run_select = Step(
        'run-select',
        template=PythonOPTemplate(
            run_select_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            slices=Slices(sub_path = True,
                input_parameter=["task_name", "trust_lvl_1", "trust_lvl_2"],
                input_artifact=["cluster_selection_index", "cluster_selection_data", "xtc_traj", "topology"],
                output_artifact=["selected_confs", "selected_cv_init", "model_devi", "selected_indices","selected_conf_tags"]
            ),
            **run_template_config,
        ),
        parameters={
            "label_config":  select_steps.inputs.parameters["label_config"],
            "trust_lvl_1": select_steps.inputs.parameters["trust_lvl_1"],
            "trust_lvl_2": select_steps.inputs.parameters["trust_lvl_2"],
            "dt": select_steps.inputs.parameters["dt"],
            "output_freq": select_steps.inputs.parameters["output_freq"],
            "slice_mode": select_steps.inputs.parameters["slice_mode"],
            "type_map": select_steps.inputs.parameters["type_map"],
            "task_name": select_steps.inputs.parameters['task_names']
        },
        artifacts={
            "cluster_selection_index": prep_select.outputs.artifacts["cluster_selection_index"],
            "cluster_selection_data": prep_select.outputs.artifacts["cluster_selection_data"],
            "models": select_steps.inputs.artifacts["models"],
            "xtc_traj": select_steps.inputs.artifacts["xtc_traj"],
            "topology": select_steps.inputs.artifacts["topology"]
        },
        key = step_keys["run_select"]+"-{{item.order}}",
        executor = run_executor,
        **run_config,
    )
    select_steps.add(run_select)

    select_steps.outputs.parameters["cluster_threshold"].value_from_parameter = prep_select.outputs.parameters["cluster_threshold"]
    select_steps.outputs.parameters["numb_cluster"].value_from_parameter = prep_select.outputs.parameters["numb_cluster"]

    select_steps.outputs.artifacts["selected_conf_tags"]._from = run_select.outputs.artifacts["selected_conf_tags"]
    select_steps.outputs.artifacts["cluster_selection_index"]._from = prep_select.outputs.artifacts["cluster_selection_index"]
    select_steps.outputs.artifacts["selected_confs"]._from = run_select.outputs.artifacts["selected_confs"]
    select_steps.outputs.artifacts["selected_cv_init"]._from = run_select.outputs.artifacts["selected_cv_init"]
    select_steps.outputs.artifacts["model_devi"]._from = run_select.outputs.artifacts["model_devi"]
    select_steps.outputs.artifacts["selected_indices"]._from = run_select.outputs.artifacts["selected_indices"]
    
    
    
    return select_steps
