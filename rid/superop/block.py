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


class FirstIterationBlock(Steps):
    def __init__(
        self,
        name: str,
        exploration_op: OP,
        select_op: OP,
        label_op: OP,
        data_op: OP,
        train_op: OP,
        train_config: Dict,
        upload_python_package
    ):

        self._input_parameters = {
            "trust_lvl_1" : InputParameter(type=float, value=2.0),
            "trust_lvl_2": InputParameter(type=float, value=3.0),
            "exploration_gmx_config" : InputParameter(type=Dict),
            "cv_config" : InputParameter(type=Dict),
            "cluster_threshold": InputParameter(type=float, value=1.0),,
            "angular_mask": InputParameter(type=Optional[Union[np.ndarray, List]]),
            "weights": InputParameter(type=Optional[Union[np.ndarray, List]]),
            "numb_cluster_upper": InputParameter(type=float),
            "numb_cluster_lower": InputParameter(type=float),
            "max_selection": InputParameter(type=int),
            "numb_cluster_threshold": InputParameter(type=float, value=30),
            "dt": InputParameter(type=float, value=0.02),
            "slice_mode": InputParameter(type=str, value="gmx"),
            "label_gmx_config": InputParameter(type=Dict),
            "kappas": InputParameter(type=List[float]),
            "angular_mask": InputParameter(type=List),
            "tail": InputParameter(type=float, value=0.9),
            "if_make_threshold": select_steps.inputs.parameters['if_make_threshold'],
            "numb_models": InputParameter(type=int, value=4),
            "train_config": InputParameter(type=Dict),
        }        
        self._input_artifacts = {
            "models" : InputArtifact(optional=True),
            "topology" : InputArtifact(),
            "confs" : InputArtifact(),
        }
        self._output_parameters = {
            "cluster_threshold": OutputParameter(type=int)
        }
        self._output_artifacts = {
            "exploration_mdrun_log": OutputArtifact(),
            "exploration_trajectory": OutputArtifact(),
            "selection_index": OutputArtifact(),
            "models": OutputArtifact()
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

        self._init_keys = ['scheduler', 'id']
        self.loop_key = 'loop'
        self.step_keys = {}
        for ii in self._init_keys:
        self.step_keys[ii] = '--'.join(['init', ii])

        self = _rid_iteration(
            self, 
            exploration_op,
            select_op,
            label_op,
            data_op,
            train_op,
            train_config,
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


def _block(
        block_steps,
        exploration_op: OP,
        select_op: OP,
        label_op: OP,
        data_op: OP,
        train_op: OP,
        train_config : Dict,
        upload_python_package : str = None,
    ):
    tarin_config = deepcopy(train_config)
    train_template_config = train_config.pop('template_config')
    train_executor = init_executor(train_config.pop('executor'))

    exploration = Step(
        'exploration',
        template=exploration_op,
        parameters={
            "trust_lvl_1" : block_steps.inputs.parameters['trust_lvl_1'],
            "trust_lvl_2": block_steps.inputs.parameters['trust_lvl_2'],
            "gmx_config" : block_steps.inputs.parameters['exploration_gmx_config'],
            "cv_config" : block_steps.inputs.parameters['cv_config']
        },
        artifacts={
            "models" : block_steps.inputs.artifacts['models'],
            "topology" : block_steps.inputs.artifacts['topology'],
            "confs" : block_steps.inputs.artifacts['confs']
        },
        key = 'exploration',
    )
    block_steps.add(exploration)

    selection = Step(
        'selection',
        template=select_op,
        parameters={
            "trust_lvl_1" : block_steps.inputs.parameters["trust_lvl_1"],
            "trust_lvl_2": block_steps.inputs.parameters["trust_lvl_2"],
            "cluster_threshold": block_steps.inputs.parameters["cluster_threshold"],
            "angular_mask": block_steps.inputs.parameters["angular_mask"],
            "weights": block_steps.inputs.parameters["weights"],
            "numb_cluster_upper": block_steps.inputs.parameters["numb_cluster_upper"],
            "numb_cluster_lower": block_steps.inputs.parameters["numb_cluster_lower"],
            "max_selection": block_steps.inputs.parameters["max_selection"],
            "numb_cluster_threshold": block_steps.inputs.parameters["numb_cluster_threshold"],
            "dt": block_steps.inputs.parameters["dt"],
            "slice_mode": block_steps.inputs.parameters["slice_mode"],
            "if_make_threshold": block_steps.inputs.parameters["if_make_threshold"],
        },
        artifacts={
            "model_list" : block_steps.inputs.artifacts["model_list"],
            "plm_out": exploration.outputs.artifacts["plm_out"],
            "xtc_traj": exploration.outputs.artifacts["xtc_traj"],
            "topology": block_steps.inputs.artifacts["topology"]
        },
        key = "selection",
    )
    block_steps.add(selection)

    label = Step(
        'label',
        template=label_op,
        parameters={
            "angular_mask": block_steps.inputs.parameters['angular_mask']
            "gmx_config": block_steps.inputs.parameters['label_gmx_config'],
            "cv_config": block_steps.inputs.parameters['cv_config'],
            "task_id": block_steps.inputs.parameters['task_id'],
            "kappas": block_steps.inputs.parameters['kappas'],
            "tail": block_steps.inputs.parameters['tail'],
        },
        artifacts={
            "topology": exploration.inputs.artifacts["plm_out"],
            "confs": selection.outputs.artifacts["selected_confs"],
            "at": selection.outputs.artifacts["selected_cv_init"]
        },
        key = "label",
        **run_config,
    )
    block_steps.add(label)

    gen_data = Step(
        'gen_data',
        template=data_op,
        parameters={},
        artifacts={
            "forces": label.outputs.artifacts["forces"],
            "centers": selection.outputs.artifacts["selected_cv_init"]
        },
        key = "gen-data",
        **run_config,
    )
    block_steps.add(gen_data)

    train = Step(
        'gen_data',
        template=PythonOPTemplate(
            train_op,
            python_packages = upload_python_package,
            **run_template_config,
        ),
        parameters={
            "task_id": "{{item}}",
            "cv_dim": exploration.outputs.parameters["cv_dim"],
            "angular_mask": block_steps.inputs.parameters["angular_mask"],
            "train_config": block_steps.inputs.parameters["train_config"],
        },
        artifacts={
            "data": gen_data.outputs.artifacts["data"],
        },
        with_param=argo_range(block_steps.inputs.parameters["numb_models"]),
        key = "gen-data",
        **run_config,
    )
    block_steps.add(train)

    block_steps.outputs.artifacts["models"]._from = train.outputs.artifacts["models"]
    block_steps.outputs.artifacts["exploration_trajectory"]._from = exploration.outputs.artifacts["trajectory"]
    block_steps.outputs.artifacts["exploration_mdrun_log"]._from = exploration.outputs.artifacts["gmx_mdrun_log"]
    block_steps.outputs.artifacts["selection_index"]._from = selection.outputs.artifacts["selected_indices"]
    block_steps.outputs.parameters["cluster_threshold"]._from = selection.outputs.artifacts["cluster_threshold"]
    
    return label_steps
