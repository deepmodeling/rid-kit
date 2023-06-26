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
    Workflow,
    Step,
    Steps,
    upload_artifact,
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


class InitBlock(Steps):
    
    r"""Initial Block SuperOP
    This SuperOP is the first iteration of the rid-kit cycle.
    """
    def __init__(
        self,
        name: str,
        exploration_op: OP,
        select_op: OP,
        label_op: OP,
        data_op: OP,
        train_op: OP,
        model_devi_op: OP,
        train_config: Dict,
        model_devi_config: Dict,
        upload_python_package = None,
        retry_times = None
    ):

        self._input_parameters = {
            "block_tag" : InputParameter(type=str, value=""),
            "walker_tags": InputParameter(type=List),
            "model_tags": InputParameter(type=List),
            "trust_lvl_1" : InputParameter(type=List[float], value=2.0),
            "trust_lvl_2": InputParameter(type=List[float], value=3.0),
            "exploration_config" : InputParameter(type=Dict),
            "cv_config" : InputParameter(type=Dict),
            "cluster_threshold": InputParameter(type=float, value=1.0),
            "angular_mask": InputParameter(type=Optional[Union[np.ndarray, List]]),
            "weights": InputParameter(type=Optional[Union[np.ndarray, List]]),
            "numb_cluster_upper": InputParameter(type=float),
            "numb_cluster_lower": InputParameter(type=float),
            "max_selection": InputParameter(type=int),
            "std_threshold": InputParameter(type=float, value=5.0),
            "dt": InputParameter(type=float, value=0.02),
            "output_freq": InputParameter(type=float, value=2500),
            "slice_mode": InputParameter(type=str, value="gmx"),
            "label_config": InputParameter(type=Dict),
            "type_map": InputParameter(type=List, value = []),
            "tail": InputParameter(type=float, value=0.9),
            "train_config": InputParameter(type=Dict)
        }        
        self._input_artifacts = {
            "models" : InputArtifact(optional=True),
            "forcefield" : InputArtifact(optional=True),
            "topology" : InputArtifact(optional=True),
            "inputfile": InputArtifact(optional=True),
            "confs" : InputArtifact(),
            "index_file": InputArtifact(optional=True),
            "data_old": InputArtifact(optional=True),
            "dp_files": InputArtifact(optional=True),
            "cv_file": InputArtifact(optional=True)
        }
        self._output_parameters = {
            "cluster_threshold": OutputParameter(type=int)
        }
        self._output_artifacts = {
            "exploration_md_log": OutputArtifact(),
            "exploration_trajectory": OutputArtifact(),
            "selection_index": OutputArtifact(),
            "models": OutputArtifact(),
            "data": OutputArtifact(),
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
        self = _first_run_block(
            self, 
            exploration_op,
            select_op,
            label_op,
            data_op,
            train_op,
            model_devi_op,
            train_config,
            model_devi_config,
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


def _first_run_block(
        block_steps,
        exploration_op: OP,
        select_op: OP,
        label_op: OP,
        data_op: OP,
        train_op: OP,
        model_devi_op: OP,
        train_config : Dict,
        model_devi_config: Dict,
        upload_python_package : str = None,
        retry_times: int = None
    ):
    exploration = Step(
        "Exploration",
        template=exploration_op,
        parameters={
            "trust_lvl_1" : block_steps.inputs.parameters['trust_lvl_1'],
            "trust_lvl_2": block_steps.inputs.parameters['trust_lvl_2'],
            "exploration_config" : block_steps.inputs.parameters['exploration_config'],
            "cv_config" : block_steps.inputs.parameters['cv_config'],
            "task_names" : block_steps.inputs.parameters['walker_tags'],
            "block_tag" : block_steps.inputs.parameters['block_tag']
        },
        artifacts={
            "models" : block_steps.inputs.artifacts['models'],
            "forcefield" : block_steps.inputs.artifacts['forcefield'],
            "topology" : block_steps.inputs.artifacts['topology'],
            "inputfile": block_steps.inputs.artifacts['inputfile'],
            "confs" : block_steps.inputs.artifacts['confs'],
            "index_file": block_steps.inputs.artifacts['index_file'],
            "dp_files": block_steps.inputs.artifacts['dp_files'],
            "cv_file": block_steps.inputs.artifacts['cv_file']
        },
        key = '{}-exploration'.format(block_steps.inputs.parameters['block_tag'])
    )
    block_steps.add(exploration)

    selection = Step(
        "Selection",
        template=select_op,
        parameters={
            "label_config": block_steps.inputs.parameters["label_config"],
            "trust_lvl_1" : block_steps.inputs.parameters["trust_lvl_1"],
            "trust_lvl_2": block_steps.inputs.parameters["trust_lvl_2"],
            "cluster_threshold": block_steps.inputs.parameters["cluster_threshold"],
            "angular_mask": block_steps.inputs.parameters["angular_mask"],
            "weights": block_steps.inputs.parameters["weights"],
            "numb_cluster_upper": block_steps.inputs.parameters["numb_cluster_upper"],
            "numb_cluster_lower": block_steps.inputs.parameters["numb_cluster_lower"],
            "max_selection": block_steps.inputs.parameters["max_selection"],
            "dt": block_steps.inputs.parameters["dt"],
            "output_freq": block_steps.inputs.parameters["output_freq"],
            "slice_mode": block_steps.inputs.parameters["slice_mode"],
            "type_map": block_steps.inputs.parameters["type_map"],
            "if_make_threshold": True,
            "task_names" : block_steps.inputs.parameters['walker_tags'],
            "block_tag" : block_steps.inputs.parameters['block_tag'],
        },
        artifacts={
            "models" : block_steps.inputs.artifacts["models"],
            "plm_out": exploration.outputs.artifacts["plm_out"],
            "xtc_traj": exploration.outputs.artifacts["trajectory"],
            "topology": block_steps.inputs.artifacts["confs"]
        },
        key = '{}-selection'.format(block_steps.inputs.parameters['block_tag']),
    )
    block_steps.add(selection)

    label = Step(
        "Label",
        template=label_op,
        parameters={
            "label_config": block_steps.inputs.parameters['label_config'],
            "cv_config": block_steps.inputs.parameters['cv_config'],
            "tail": block_steps.inputs.parameters['tail'],
            "block_tag" : block_steps.inputs.parameters['block_tag'],
            "std_threshold": block_steps.inputs.parameters["std_threshold"]
        },
        artifacts={
            "topology": block_steps.inputs.artifacts["topology"],
            "forcefield" : block_steps.inputs.artifacts['forcefield'],
            "confs": selection.outputs.artifacts["selected_confs"],
            "at": selection.outputs.artifacts["selected_cv_init"],
            "index_file": block_steps.inputs.artifacts['index_file'],
            "inputfile": block_steps.inputs.artifacts['inputfile'],
            "dp_files": block_steps.inputs.artifacts['dp_files'],
            "cv_file": block_steps.inputs.artifacts['cv_file'],
            "conf_tags" : selection.outputs.artifacts['selected_conf_tags']
        },
        key = '{}-label'.format(block_steps.inputs.parameters['block_tag'])
    )
    block_steps.add(label)

    gen_data = Step(
        'GenData',
        template=data_op,
        parameters={"block_tag" : block_steps.inputs.parameters['block_tag']},
        artifacts={
            "cv_forces": label.outputs.artifacts["cv_forces"],
            "data_old": block_steps.inputs.artifacts['data_old']
        },
        key = '{}-gen-data'.format(block_steps.inputs.parameters['block_tag']),
    )
    block_steps.add(gen_data)

    train_config = deepcopy(train_config)
    train_template_config = train_config.pop('template_config')
    train_executor = init_executor(train_config.pop('executor'))
    train = Step(
        "train",
        template=PythonOPTemplate(
            train_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            slices=Slices("{{item}}",
                input_parameter=["model_tag"],
                output_artifact=["model","train_fig"]),
            **train_template_config,
        ),
        parameters={
            "model_tag": block_steps.inputs.parameters["model_tags"],
            "angular_mask": block_steps.inputs.parameters["angular_mask"],
            "train_config": block_steps.inputs.parameters["train_config"],
        },
        artifacts={
            "data": gen_data.outputs.artifacts["data"],
        },
        executor = train_executor,
        with_param=argo_range(argo_len(block_steps.inputs.parameters["model_tags"])),
        key = "{}-train".format(block_steps.inputs.parameters["block_tag"])+"-{{item}}",
        **train_config,
    )
    block_steps.add(train)
    
    model_devi_config = deepcopy(model_devi_config)
    model_devi_template_config = model_devi_config.pop('template_config')
    model_devi_executor = init_executor(model_devi_config.pop('executor'))
    deviation = Step(
        "ModelDeviation",
        template=PythonOPTemplate(
            model_devi_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            slices=Slices("{{item}}",
                input_parameter=["trust_lvl_1","task_name"],
                input_artifact=["plm_out","selected_indices"],
                output_artifact=["model_devi","model_devi_fig"]),
            **model_devi_template_config,
        ),
        parameters={
            "trust_lvl_1": block_steps.inputs.parameters["trust_lvl_1"],
            "block_tag" : block_steps.inputs.parameters['block_tag'],
            "task_name": block_steps.inputs.parameters['walker_tags'],
        },
        artifacts={
            "models" : train.outputs.artifacts["model"],
            "plm_out": exploration.outputs.artifacts["plm_out"],
            "selected_indices": selection.outputs.artifacts["selected_indices"]
        },
        executor = model_devi_executor,
        with_param=argo_range(argo_len(block_steps.inputs.parameters["walker_tags"])),
        key = '{}-model-devi'.format(block_steps.inputs.parameters['block_tag'])+"-{{item}}",
        **model_devi_config
    )
    block_steps.add(deviation)

    block_steps.outputs.artifacts["models"]._from = train.outputs.artifacts["model"]
    block_steps.outputs.artifacts["data"]._from = gen_data.outputs.artifacts["data"]
    block_steps.outputs.artifacts["conf_outs"]._from = exploration.outputs.artifacts["conf_outs"]
    block_steps.outputs.artifacts["exploration_trajectory"]._from = exploration.outputs.artifacts["trajectory"]
    block_steps.outputs.artifacts["exploration_md_log"]._from = exploration.outputs.artifacts["md_log"]
    block_steps.outputs.artifacts["selection_index"]._from = selection.outputs.artifacts["selected_indices"]
    block_steps.outputs.parameters["cluster_threshold"].value_from_parameter = selection.outputs.parameters["cluster_threshold"]
    
    return block_steps


class IterBlock(Steps):
    
    r"""Iterative Block SuperOP.
    This SuperOP is the iterations after the inital iteration of the rid-kit cycle.
    """
    def __init__(
        self,
        name: str,
        exploration_op: OP,
        select_op: OP,
        label_op: OP,
        data_op: OP,
        adjust_lvl_op: OP,
        train_op: OP,  
        model_devi_op: OP,
        adjust_lvl_config: Dict,
        train_config: Dict,
        model_devi_config: Dict,
        upload_python_package = None,
        retry_times = None
    ):

        self._input_parameters = {
            "block_tag" : InputParameter(type=str, value=""),
            "walker_tags": InputParameter(type=List),
            "model_tags": InputParameter(type=List),
            "trust_lvl_1" : InputParameter(type=List[float]),
            "trust_lvl_2": InputParameter(type=List[float]),
            "init_trust_lvl_1" : InputParameter(type=List[float]),
            "init_trust_lvl_2": InputParameter(type=List[float]),
            "exploration_config" : InputParameter(type=Dict),
            "cv_config" : InputParameter(type=Dict),
            "cluster_threshold": InputParameter(type=float, value=1.0),
            "angular_mask": InputParameter(type=Optional[Union[np.ndarray, List]]),
            "weights": InputParameter(type=Optional[Union[np.ndarray, List]]),
            "max_selection": InputParameter(type=int),
            "numb_cluster_threshold": InputParameter(type=float, value=30),
            "std_threshold": InputParameter(type=float, value=5.0),
            "dt": InputParameter(type=float, value=0.02),
            "output_freq": InputParameter(type=float, value=2500),
            "slice_mode": InputParameter(type=str, value="gmx"),
            "label_config": InputParameter(type=Dict),
            "tail": InputParameter(type=float, value=0.9),
            "train_config": InputParameter(type=Dict),
            "type_map": InputParameter(type=List, value=[]),
            "adjust_amplifier": InputParameter(type=float, value=1.5),
            "max_level_multiple": InputParameter(type=float, value=8.0),
        }        
        self._input_artifacts = {
            "models" : InputArtifact(optional=True),
            "forcefield" : InputArtifact(optional=True),
            "topology" : InputArtifact(optional=True),
            "inputfile": InputArtifact(optional=True),
            "confs" : InputArtifact(),
            "data_old": InputArtifact(),
            "index_file": InputArtifact(optional=True),
            "dp_files": InputArtifact(optional=True),
            "cv_file": InputArtifact(optional=True)
        }
        self._output_parameters = {
            "cluster_threshold": OutputParameter(type=int),
            "adjust_trust_lvl_1": OutputParameter(type=int),
            "adjust_trust_lvl_2": OutputParameter(type=int),
        }
        self._output_artifacts = {
            "exploration_md_log": OutputArtifact(),
            "exploration_trajectory": OutputArtifact(),
            "selection_index": OutputArtifact(),
            "models": OutputArtifact(),
            "data": OutputArtifact(),
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

        self = _iter_block(
            self, 
            exploration_op,
            select_op,
            label_op,
            data_op,
            adjust_lvl_op,
            train_op,
            model_devi_op,
            adjust_lvl_config,
            train_config,
            model_devi_config,
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


def _iter_block(
        block_steps,
        exploration_op: OP,
        select_op: OP,
        label_op: OP,
        data_op: OP,
        adjust_lvl_op: OP,
        train_op: OP,
        model_devi_op: OP,
        adjust_lvl_config : Dict,
        train_config : Dict,
        model_devi_config: Dict,
        upload_python_package : str = None,
        retry_times: int = None
    ):

    exploration = Step(
        "Exploration",
        template=exploration_op,
        parameters={
            "trust_lvl_1" : block_steps.inputs.parameters['trust_lvl_1'],
            "trust_lvl_2": block_steps.inputs.parameters['trust_lvl_2'],
            "exploration_config" : block_steps.inputs.parameters['exploration_config'],
            "cv_config" : block_steps.inputs.parameters['cv_config'],
            "task_names" : block_steps.inputs.parameters['walker_tags'],
            "block_tag" : block_steps.inputs.parameters['block_tag'],
        },
        artifacts={
            "models" : block_steps.inputs.artifacts['models'],
            "forcefield" : block_steps.inputs.artifacts['forcefield'],
            "topology" : block_steps.inputs.artifacts['topology'],
            "inputfile": block_steps.inputs.artifacts['inputfile'],
            "confs" : block_steps.inputs.artifacts['confs'],
            "index_file": block_steps.inputs.artifacts['index_file'],
            "dp_files": block_steps.inputs.artifacts['dp_files'],
            "cv_file": block_steps.inputs.artifacts['cv_file']
        },
        key = '{}-exploration'.format(block_steps.inputs.parameters['block_tag'])
    )
    block_steps.add(exploration)

    selection = Step(
        "Selection",
        template=select_op,
        parameters={
            "label_config": block_steps.inputs.parameters["label_config"],
            "trust_lvl_1" : block_steps.inputs.parameters["init_trust_lvl_1"],
            "trust_lvl_2": block_steps.inputs.parameters["init_trust_lvl_2"],
            "cluster_threshold": block_steps.inputs.parameters["cluster_threshold"],
            "angular_mask": block_steps.inputs.parameters["angular_mask"],
            "weights": block_steps.inputs.parameters["weights"],
            "max_selection": block_steps.inputs.parameters["max_selection"],
            "dt": block_steps.inputs.parameters["dt"],
            "output_freq": block_steps.inputs.parameters["output_freq"],
            "slice_mode": block_steps.inputs.parameters["slice_mode"],
            "type_map": block_steps.inputs.parameters["type_map"],
            "if_make_threshold": False,
            "task_names" : block_steps.inputs.parameters['walker_tags'],
            "block_tag" : block_steps.inputs.parameters['block_tag'],
        },
        artifacts={
            "models" : block_steps.inputs.artifacts["models"],
            "plm_out": exploration.outputs.artifacts["plm_out"],
            "xtc_traj": exploration.outputs.artifacts["trajectory"],
            "topology": block_steps.inputs.artifacts["confs"]
        },
        key = '{}-selection'.format(block_steps.inputs.parameters['block_tag']),
    )
    block_steps.add(selection)

    adjust_lvl_config = deepcopy(adjust_lvl_config)
    adjust_lvl_template_config = adjust_lvl_config.pop('template_config')
    adjust_lvl_executor = init_executor(adjust_lvl_config.pop('executor'))

    adjust_lvl = Step(
        "adjust-level",
        template=PythonOPTemplate(
            adjust_lvl_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            slices=Slices("{{item}}",
                input_parameter=["trust_lvl_1", "trust_lvl_2", "numb_cluster", "init_trust_lvl_1", "init_trust_lvl_2"],
                output_parameter=["adjust_trust_lvl_1", "adjust_trust_lvl_2"]
                ),
            **adjust_lvl_template_config,
        ),
        parameters={
            "trust_lvl_1": block_steps.inputs.parameters["trust_lvl_1"],
            "trust_lvl_2": block_steps.inputs.parameters["trust_lvl_2"],
            "init_trust_lvl_1": block_steps.inputs.parameters["init_trust_lvl_1"],
            "init_trust_lvl_2": block_steps.inputs.parameters["init_trust_lvl_2"],
            "numb_cluster": selection.outputs.parameters["numb_cluster"],
            "numb_cluster_threshold": block_steps.inputs.parameters["numb_cluster_threshold"],
            "adjust_amplifier": block_steps.inputs.parameters["adjust_amplifier"], 
            "max_level_multiple": block_steps.inputs.parameters["max_level_multiple"]
        },
        artifacts={},
        with_param=argo_range(argo_len(block_steps.inputs.parameters["trust_lvl_1"])),
        executor = adjust_lvl_executor,
        key = '{}-adjust-level'.format(block_steps.inputs.parameters['block_tag'])+"-{{item}}",
        **adjust_lvl_config,
    )
    block_steps.add(adjust_lvl)

    label = Step(
        "Label",
        template=label_op,
        parameters={
            "label_config": block_steps.inputs.parameters['label_config'],
            "cv_config": block_steps.inputs.parameters['cv_config'],
            "tail": block_steps.inputs.parameters['tail'],
            "block_tag" : block_steps.inputs.parameters['block_tag'],
            "std_threshold": block_steps.inputs.parameters["std_threshold"]
        },
        artifacts={
            "topology": block_steps.inputs.artifacts["topology"],
            "forcefield" : block_steps.inputs.artifacts['forcefield'],
            "confs": selection.outputs.artifacts["selected_confs"],
            "inputfile": block_steps.inputs.artifacts['inputfile'],
            "at": selection.outputs.artifacts["selected_cv_init"],
            "index_file": block_steps.inputs.artifacts['index_file'],
            "dp_files": block_steps.inputs.artifacts['dp_files'],
            "cv_file": block_steps.inputs.artifacts['cv_file'],
            "conf_tags" : selection.outputs.artifacts['selected_conf_tags']
        },
        key = '{}-label'.format(block_steps.inputs.parameters['block_tag'])
    )
    block_steps.add(label)

    gen_data = Step(
        'GenData',
        template=data_op,
        parameters={"block_tag" : block_steps.inputs.parameters['block_tag']},
        artifacts={
            "cv_forces": label.outputs.artifacts["cv_forces"],
            "data_old": block_steps.inputs.artifacts['data_old']
        },
        key = '{}-gen-data'.format(block_steps.inputs.parameters['block_tag']),
    )
    block_steps.add(gen_data)

    train_config = deepcopy(train_config)
    train_template_config = train_config.pop('template_config')
    train_executor = init_executor(train_config.pop('executor'))
    train = Step(
        "train",
        template=PythonOPTemplate(
            train_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            slices=Slices("{{item}}",
                input_parameter=["model_tag"],
                output_artifact=["model","train_fig"]),
            **train_template_config,
        ),
        parameters={
            "model_tag": block_steps.inputs.parameters["model_tags"],
            "angular_mask": block_steps.inputs.parameters["angular_mask"],
            "train_config": block_steps.inputs.parameters["train_config"],
        },
        artifacts={
            "data": gen_data.outputs.artifacts["data"],
        },
        executor = train_executor,
        with_param=argo_range(argo_len(block_steps.inputs.parameters["model_tags"])),
        key = "{}-train".format(block_steps.inputs.parameters["block_tag"])+"-{{item}}",
        **train_config,
    )
    block_steps.add(train)
    
    model_devi_config = deepcopy(model_devi_config)
    model_devi_template_config = model_devi_config.pop('template_config')
    model_devi_executor = init_executor(model_devi_config.pop('executor'))
    deviation = Step(
        "ModelDeviation",
        template=PythonOPTemplate(
            model_devi_op,
            python_packages = upload_python_package,
            retry_on_transient_error = retry_times,
            slices=Slices("{{item}}",
                input_parameter=["trust_lvl_1","task_name"],
                input_artifact=["plm_out","selected_indices"],
                output_artifact=["model_devi","model_devi_fig"]),
            **model_devi_template_config,
        ),
        parameters={
            "trust_lvl_1": block_steps.inputs.parameters["trust_lvl_1"],
            "block_tag" : block_steps.inputs.parameters['block_tag'],
            "task_name": block_steps.inputs.parameters['walker_tags'],
        },
        artifacts={
            "models" : train.outputs.artifacts["model"],
            "plm_out": exploration.outputs.artifacts["plm_out"],
            "selected_indices": selection.outputs.artifacts["selected_indices"]
        },
        executor = model_devi_executor,
        with_param=argo_range(argo_len(block_steps.inputs.parameters["walker_tags"])),
        key = '{}-model-devi'.format(block_steps.inputs.parameters['block_tag'])+"-{{item}}",
        **model_devi_config
    )
    block_steps.add(deviation)

    block_steps.outputs.artifacts["models"]._from = train.outputs.artifacts["model"]
    block_steps.outputs.artifacts["data"]._from = gen_data.outputs.artifacts["data"]
    block_steps.outputs.artifacts["conf_outs"]._from = exploration.outputs.artifacts["conf_outs"]
    block_steps.outputs.artifacts["exploration_trajectory"]._from = exploration.outputs.artifacts["trajectory"]
    block_steps.outputs.artifacts["exploration_md_log"]._from = exploration.outputs.artifacts["md_log"]
    block_steps.outputs.artifacts["selection_index"]._from = selection.outputs.artifacts["selected_indices"]
    block_steps.outputs.parameters["cluster_threshold"].value_from_parameter = selection.outputs.parameters["cluster_threshold"]
    block_steps.outputs.parameters["adjust_trust_lvl_1"].value_from_parameter = adjust_lvl.outputs.parameters["adjust_trust_lvl_1"]
    block_steps.outputs.parameters["adjust_trust_lvl_2"].value_from_parameter = adjust_lvl.outputs.parameters["adjust_trust_lvl_2"]
    
    return block_steps
