from dflow import (
    InputParameter,
    Inputs,
    InputArtifact,
    Outputs,
    OutputArtifact,
    Step,
    Steps,
    if_expression,
)
from dflow.python import(
    PythonOPTemplate,
    OP
)
from typing import List, Optional, Dict, Union
import numpy as np
from rid.utils import init_executor
from rid.op.prep_rid import PrepRiD
from rid.op.recorder import Recorder
from copy import deepcopy


class ReinforcedDynamicsLoop(Steps):
    def __init__(
            self,
            name : str,
            block_op : Steps,
            step_config : dict,
            upload_python_package : str = None,
    ):
        self._input_parameters={
            "numb_iters" : InputParameter(type=int),
            "last_iteration" : InputParameter(type=int),
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
            "type_map": InputParameter(type=List, value = []),
            "adjust_amplifier": InputParameter(type=float, value=1.5),
            "max_level_multiple": InputParameter(type=float, value=8.0),
        }
        self._input_artifacts={
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
        self._output_parameters={
        }
        self._output_artifacts={
            "exploration_trajectory": OutputArtifact(),
            "models": OutputArtifact(),
            "data": OutputArtifact(),
            "conf_outs": OutputArtifact()
        }
        
        super().__init__(
            name = name,
            inputs = Inputs(
                parameters=self._input_parameters,
                artifacts=self._input_artifacts,
            ),
            outputs=Outputs(
                parameters=self._output_parameters,
                artifacts=self._output_artifacts,
            ),
        )

        _step_keys = ['block', 'recorder']
        step_keys = {}
        for ii in _step_keys:
            step_keys[ii] = '-'.join(["%s"%self.inputs.parameters["block_tag"], ii])
        
        self = _loop(
            self,
            step_keys,
            name,
            block_op,
            step_config = step_config,
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


def _loop (
        steps, 
        step_keys,
        name : str,
        block_op : OP,
        step_config : dict,
        upload_python_package : str = None,
):    
    step_config = deepcopy(step_config)
    step_template_config = step_config.pop('template_config')
    step_executor = init_executor(step_config.pop('executor'))

    recorder_step = Step(
        name = name + '-recorder',
        template=PythonOPTemplate(
            Recorder,
            python_packages = upload_python_package,
            **step_template_config,
        ),
        parameters={
            "iteration": steps.inputs.parameters['last_iteration'],
        },
        artifacts={},
        key = step_keys['recorder'],
        executor = step_executor,
        **step_config,
    )
    steps.add(recorder_step)

    block_step = Step(
        name = name + '-block',
        template = block_op,
        parameters={
            "block_tag": recorder_step.outputs.parameters["block_tag"],
            "walker_tags": steps.inputs.parameters["walker_tags"],
            "model_tags": steps.inputs.parameters["model_tags"],
            "exploration_config": steps.inputs.parameters["exploration_config"],
            "cv_config": steps.inputs.parameters["cv_config"],
            "trust_lvl_1" : steps.inputs.parameters["trust_lvl_1"],
            "trust_lvl_2": steps.inputs.parameters["trust_lvl_2"],
            "init_trust_lvl_1": steps.inputs.parameters["init_trust_lvl_1"],
            "init_trust_lvl_2": steps.inputs.parameters["init_trust_lvl_2"],
            "cluster_threshold": steps.inputs.parameters["cluster_threshold"],
            "angular_mask": steps.inputs.parameters["angular_mask"],
            "weights": steps.inputs.parameters["weights"],
            "max_selection": steps.inputs.parameters["max_selection"],
            "numb_cluster_threshold": steps.inputs.parameters["numb_cluster_threshold"],
            "std_threshold": steps.inputs.parameters["std_threshold"],
            "dt": steps.inputs.parameters["dt"],
            "output_freq": steps.inputs.parameters["output_freq"],
            "slice_mode": steps.inputs.parameters["slice_mode"],
            "type_map": steps.inputs.parameters["type_map"],
            "label_config": steps.inputs.parameters["label_config"],
            "train_config": steps.inputs.parameters["train_config"]
        },
        artifacts={
            "models": steps.inputs.artifacts["models"],
            "forcefield" : steps.inputs.artifacts['forcefield'],
            "topology": steps.inputs.artifacts["topology"],
            "inputfile": steps.inputs.artifacts["inputfile"],
            "confs": steps.inputs.artifacts["confs"],
            "data_old": steps.inputs.artifacts["data_old"],
            "index_file": steps.inputs.artifacts["index_file"],
            "dp_files": steps.inputs.artifacts["dp_files"],
            "cv_file": steps.inputs.artifacts["cv_file"]
        },
        key = step_keys['block'],
    )
    steps.add(block_step)

    next_step = Step(
        name = name+'-next',
        template = steps,
        parameters={
            "numb_iters": steps.inputs.parameters["numb_iters"],
            "last_iteration": recorder_step.outputs.parameters["next_iteration"],
            "block_tag": recorder_step.outputs.parameters["block_tag"],
            "walker_tags": steps.inputs.parameters["walker_tags"],
            "model_tags": steps.inputs.parameters["model_tags"],
            "exploration_config": steps.inputs.parameters["exploration_config"],
            "cv_config": steps.inputs.parameters["cv_config"],
            "trust_lvl_1": block_step.outputs.parameters["adjust_trust_lvl_1"],
            "trust_lvl_2": block_step.outputs.parameters["adjust_trust_lvl_2"],
            "init_trust_lvl_1": steps.inputs.parameters["init_trust_lvl_1"],
            "init_trust_lvl_2": steps.inputs.parameters["init_trust_lvl_2"],
            "cluster_threshold":  block_step.outputs.parameters["cluster_threshold"],
            "angular_mask": steps.inputs.parameters["angular_mask"],
            "weights": steps.inputs.parameters["weights"],
            "max_selection": steps.inputs.parameters["max_selection"],
            "numb_cluster_threshold": steps.inputs.parameters["numb_cluster_threshold"],
            "std_threshold": steps.inputs.parameters["std_threshold"],
            "dt": steps.inputs.parameters["dt"],
            "output_freq": steps.inputs.parameters["output_freq"],
            "slice_mode": steps.inputs.parameters["slice_mode"],
            "type_map": steps.inputs.parameters["type_map"],
            "label_config": steps.inputs.parameters["label_config"],
            "train_config": steps.inputs.parameters["train_config"]
        },
        artifacts={
            "models":  block_step.outputs.artifacts["models"],
            "forcefield" : steps.inputs.artifacts['forcefield'],
            "topology": steps.inputs.artifacts["topology"],
            "inputfile": steps.inputs.artifacts["inputfile"],
            "confs": block_step.outputs.artifacts["conf_outs"],
            "data_old": block_step.outputs.artifacts["data"],
            "index_file": steps.inputs.artifacts["index_file"],
            "dp_files": steps.inputs.artifacts["dp_files"],
            "cv_file": steps.inputs.artifacts["cv_file"]
        },
        when = "%s < %s" % (recorder_step.outputs.parameters['next_iteration'], steps.inputs.parameters["numb_iters"]),
    )
    steps.add(next_step)    

    steps.outputs.artifacts['exploration_trajectory'].from_expression = \
        if_expression(
            _if = (recorder_step.outputs.parameters['next_iteration'] >= steps.inputs.parameters["numb_iters"]),
            _then = block_step.outputs.artifacts['exploration_trajectory'],
            _else = next_step.outputs.artifacts['exploration_trajectory'],
        )
    steps.outputs.artifacts['conf_outs'].from_expression = \
        if_expression(
            _if = (recorder_step.outputs.parameters['next_iteration'] >= steps.inputs.parameters["numb_iters"]),
            _then = block_step.outputs.artifacts['conf_outs'],
            _else = next_step.outputs.artifacts['conf_outs'],
        )
    steps.outputs.artifacts['models'].from_expression = \
        if_expression(
            _if = (recorder_step.outputs.parameters['next_iteration'] >= steps.inputs.parameters["numb_iters"]),
            _then = block_step.outputs.artifacts['models'],
            _else = next_step.outputs.artifacts['models'],
        )
    steps.outputs.artifacts['data'].from_expression = \
        if_expression(
            _if = (recorder_step.outputs.parameters['next_iteration'] >= steps.inputs.parameters["numb_iters"]),
            _then = block_step.outputs.artifacts['data'],
            _else = next_step.outputs.artifacts['data'],
        )

    return steps


class ReinforcedDynamics(Steps):
    def __init__(
            self,
            name : str,
            init_block_op : Steps,
            block_op : Steps,
            step_config : dict,
            upload_python_package : str = None,
    ):
        
        self._input_parameters={}
        self._input_artifacts={
            "models": InputArtifact(optional=True),
            "forcefield": InputArtifact(optional=True),
            "topology": InputArtifact(optional=True),
            "inputfile": InputArtifact(optional=True),
            "confs": InputArtifact(),
            "rid_config": InputArtifact(),
            "index_file": InputArtifact(optional=True),
            "data_file": InputArtifact(optional=True),
            "dp_files": InputArtifact(optional=True),
            "cv_file": InputArtifact(optional=True)
        }
        self._output_parameters={
        }
        self._output_artifacts={
            "exploration_trajectory": OutputArtifact(),
            "models": OutputArtifact(),
            "data": OutputArtifact(),
            "conf_outs": OutputArtifact()
        }        
        
        super().__init__(
            name = name,
            inputs = Inputs(
                parameters=self._input_parameters,
                artifacts=self._input_artifacts,
            ),
            outputs=Outputs(
                parameters=self._output_parameters,
                artifacts=self._output_artifacts,
            ),
        )

        _init_keys = ['recorder', 'block']
        step_keys = {}
        for ii in _init_keys:
            step_keys[ii] = '--'.join(['init', ii])

        self = _rid(
            self,
            step_keys,
            name, 
            init_block_op,
            block_op,
            step_config = step_config,
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
    def init_keys(self):
        return self._init_keys

    @property
    def loop_keys(self):
        return [self.loop_key] + self.loop.keys


def _rid(
        steps, 
        step_keys,
        name,
        init_block_op,
        block_op,
        step_config : dict,
        upload_python_package : Optional[str] = None
    ):  

    _step_config = deepcopy(step_config)
    step_template_config = _step_config.pop('template_config')
    step_executor = init_executor(_step_config.pop('executor'))

    prep_rid = Step(
        name = 'prepare-rid',
        template=PythonOPTemplate(
            PrepRiD,
            python_packages = upload_python_package,
            **step_template_config,
        ),
        parameters={},
        artifacts={
            "confs": steps.inputs.artifacts['confs'],
            "rid_config" : steps.inputs.artifacts['rid_config'],
        },
        key = 'prepare-rid',
        executor = step_executor,
        **_step_config,
    )
    steps.add(prep_rid)

    recorder_step = Step(
        name = name + '-recorder',
        template=PythonOPTemplate(
            Recorder,
            python_packages = upload_python_package,
            **step_template_config,
        ),
        parameters={
            "iteration": None,
        },
        artifacts={},
        key = "init-recorder",
        executor = step_executor,
        **_step_config,
    )
    steps.add(recorder_step)

    init_block = Step(
        name = name + '-block',
        template = init_block_op,
        parameters={
            "block_tag": recorder_step.outputs.parameters["block_tag"],
            "walker_tags": prep_rid.outputs.parameters["walker_tags"],
            "model_tags": prep_rid.outputs.parameters["model_tags"],
            "exploration_config": prep_rid.outputs.parameters["exploration_config"],
            "cv_config": prep_rid.outputs.parameters["cv_config"],
            "trust_lvl_1" : prep_rid.outputs.parameters["trust_lvl_1"],
            "trust_lvl_2": prep_rid.outputs.parameters["trust_lvl_2"],
            "numb_cluster_upper": prep_rid.outputs.parameters['numb_cluster_upper'],
            "numb_cluster_lower": prep_rid.outputs.parameters['numb_cluster_lower'],
            "cluster_threshold": prep_rid.outputs.parameters["cluster_threshold"],
            "angular_mask": prep_rid.outputs.parameters["angular_mask"],
            "weights": prep_rid.outputs.parameters["weights"],
            "max_selection": prep_rid.outputs.parameters["max_selection"],
            "std_threshold": prep_rid.outputs.parameters["std_threshold"],
            "dt": prep_rid.outputs.parameters["dt"],
            "output_freq": prep_rid.outputs.parameters["output_freq"],
            "slice_mode": prep_rid.outputs.parameters["slice_mode"],
            "type_map": prep_rid.outputs.parameters["type_map"],
            "label_config": prep_rid.outputs.parameters["label_config"],
            "train_config": prep_rid.outputs.parameters["train_config"]
        },
        artifacts={
            "models": steps.inputs.artifacts["models"],
            "forcefield" : steps.inputs.artifacts['forcefield'],
            "topology": steps.inputs.artifacts["topology"],
            "inputfile": steps.inputs.artifacts["inputfile"],
            "confs": prep_rid.outputs.artifacts["confs"],
            "index_file": steps.inputs.artifacts["index_file"],
            "data_old": steps.inputs.artifacts["data_file"],
            "dp_files": steps.inputs.artifacts["dp_files"],
            "cv_file": steps.inputs.artifacts["cv_file"]
        },
        key = "%s-init-block"%recorder_step.outputs.parameters["block_tag"]
    )
    steps.add(init_block)

    loop_step = Step(
        name = name + '-loop',
        template = ReinforcedDynamicsLoop(
            name = "RiD-Loop",
            block_op = block_op,
            step_config = step_config,
            upload_python_package = upload_python_package
        ),
        parameters={
            "numb_iters": prep_rid.outputs.parameters["numb_iters"],
            "last_iteration": recorder_step.outputs.parameters["next_iteration"],
            "block_tag": recorder_step.outputs.parameters["block_tag"],
            "walker_tags": prep_rid.outputs.parameters["walker_tags"],
            "model_tags": prep_rid.outputs.parameters["model_tags"],
            "exploration_config": prep_rid.outputs.parameters["exploration_config"],
            "cv_config": prep_rid.outputs.parameters["cv_config"],
            "trust_lvl_1" : prep_rid.outputs.parameters["trust_lvl_1"],
            "trust_lvl_2": prep_rid.outputs.parameters["trust_lvl_2"],
            "init_trust_lvl_1": prep_rid.outputs.parameters["trust_lvl_1"],
            "init_trust_lvl_2": prep_rid.outputs.parameters["trust_lvl_2"],
            "cluster_threshold":  init_block.outputs.parameters["cluster_threshold"],
            "angular_mask": prep_rid.outputs.parameters["angular_mask"],
            "weights": prep_rid.outputs.parameters["weights"],
            "max_selection": prep_rid.outputs.parameters["max_selection"],
            "numb_cluster_threshold": prep_rid.outputs.parameters["numb_cluster_threshold"],
            "std_threshold": prep_rid.outputs.parameters["std_threshold"],
            "dt": prep_rid.outputs.parameters["dt"],
            "output_freq": prep_rid.outputs.parameters["output_freq"],
            "slice_mode": prep_rid.outputs.parameters["slice_mode"],
            "type_map": prep_rid.outputs.parameters["type_map"],
            "label_config": prep_rid.outputs.parameters["label_config"],
            "train_config": prep_rid.outputs.parameters["train_config"]
        },
        artifacts={
            "models":  init_block.outputs.artifacts["models"],
            "forcefield" : steps.inputs.artifacts['forcefield'],
            "topology": steps.inputs.artifacts["topology"],
            "inputfile": steps.inputs.artifacts["inputfile"],
            "confs": init_block.outputs.artifacts["conf_outs"],
            "data_old": init_block.outputs.artifacts["data"],
            "index_file": steps.inputs.artifacts["index_file"],
            "dp_files": steps.inputs.artifacts["dp_files"],
            "cv_file": steps.inputs.artifacts["cv_file"]
        },
        when = "%s < %s" % (recorder_step.outputs.parameters['next_iteration'], prep_rid.outputs.parameters["numb_iters"]),
        key = "rid-loop",
    )
    steps.add(loop_step)

    steps.outputs.artifacts['exploration_trajectory'].from_expression = \
        if_expression(
            _if = (recorder_step.outputs.parameters['next_iteration'] >= prep_rid.outputs.parameters["numb_iters"]),
            _then = init_block.outputs.artifacts['exploration_trajectory'],
            _else = loop_step.outputs.artifacts['exploration_trajectory'],
        )
    steps.outputs.artifacts['conf_outs'].from_expression = \
        if_expression(
            _if = (recorder_step.outputs.parameters['next_iteration'] >= prep_rid.outputs.parameters["numb_iters"]),
            _then = init_block.outputs.artifacts['conf_outs'],
            _else = loop_step.outputs.artifacts['conf_outs'],
        )
    steps.outputs.artifacts['models'].from_expression = \
        if_expression(
            _if = (recorder_step.outputs.parameters['next_iteration'] >= prep_rid.outputs.parameters["numb_iters"]),
            _then = init_block.outputs.artifacts['models'],
            _else = loop_step.outputs.artifacts['models'],
        )
    steps.outputs.artifacts['data'].from_expression = \
        if_expression(
            _if = (recorder_step.outputs.parameters['next_iteration'] >= prep_rid.outputs.parameters["numb_iters"]),
            _then = init_block.outputs.artifacts['data'],
            _else = loop_step.outputs.artifacts['data'],
        )
    
    return steps