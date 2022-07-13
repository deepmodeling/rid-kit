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
    if_expression,
)
from dflow.python import(
    PythonOPTemplate,
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Slices,
    BigParameter,
)
import pickle, jsonpickle, os
from typing import (
    List
)
from pathlib import Path
from dpgen2.exploration.scheduler import ExplorationScheduler
from dpgen2.exploration.report import ExplorationReport
from dpgen2.exploration.task import ExplorationTaskGroup
from dpgen2.exploration.selector import ConfSelector
from dpgen2.superop.block import ConcurrentLearningBlock
from dpgen2.utils import (
    load_object_from_file,
    dump_object_to_file,
)
from dpgen2.utils.step_config import normalize as normalize_step_dict
from dpgen2.utils.step_config import init_executor

from copy import deepcopy


class SchedulerWrapper(OP):

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            "exploration_scheduler" : Artifact(Path),
            "exploration_report": BigParameter(ExplorationReport),
            "trajs": Artifact(List[Path]),
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            "exploration_scheduler" : Artifact(Path),
            "converged" : bool,
            "lmp_task_grp" : Artifact(Path),
            "conf_selector" : ConfSelector,
        })

    @OP.exec_sign_check
    def execute(
            self,
            ip : OPIO,
    ) -> OPIO:
        scheduler_in = ip['exploration_scheduler']
        report = ip['exploration_report']
        trajs = ip['trajs']
        lmp_task_grp_file = Path('lmp_task_grp.dat')
        scheduler_file = Path('scheduler.dat')

        scheduler = load_object_from_file(scheduler_in)

        conv, lmp_task_grp, selector = scheduler.plan_next_iteration(report, trajs)

        dump_object_to_file(lmp_task_grp, lmp_task_grp_file)
        dump_object_to_file(scheduler, scheduler_file)

        return OPIO({
            "exploration_scheduler" : scheduler_file,
            "converged" : conv,
            "conf_selector" : selector,
            "lmp_task_grp" : lmp_task_grp_file,
        })


class MakeBlockId(OP):
    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            "exploration_scheduler" : Artifact(Path),
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            "block_id" : str,
        })

    @OP.exec_sign_check
    def execute(
            self,
            ip : OPIO,
    ) -> OPIO:
        scheduler_in = ip['exploration_scheduler']        
        scheduler = load_object_from_file(scheduler_in)

        stage = scheduler.get_stage()
        iteration = scheduler.get_iteration()

        return OPIO({
            "block_id" : f'iter-{iteration:06d}',
        })


class ConcurrentLearningLoop(Steps):
    def __init__(
            self,
            name : str,
            block_op : Steps,
            step_config : dict = normalize_step_dict({}),
            upload_python_package : str = None,
    ):
        self._input_parameters={
            "block_id" : InputParameter(),
            "type_map" : InputParameter(),
            "numb_models": InputParameter(type=int),
            "template_script" : InputParameter(),
            "train_config" : InputParameter(),
            "lmp_config" : InputParameter(),
            "conf_selector" : InputParameter(),
            "fp_config" : InputParameter(),
        }
        self._input_artifacts={
            "exploration_scheduler" : InputArtifact(),
            "init_models" : InputArtifact(optional=True),
            "init_data" : InputArtifact(),
            "iter_data" : InputArtifact(),
            "lmp_task_grp" : InputArtifact(),
            "fp_inputs" : InputArtifact(),
        }
        self._output_parameters={
        }
        self._output_artifacts={
            "exploration_scheduler": OutputArtifact(),
            "models": OutputArtifact(),
            "iter_data" : OutputArtifact(),
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

        self._my_keys = ['block', 'scheduler', 'id']
        self._keys = \
            self._my_keys[:1] + \
            block_op.keys + \
            self._my_keys[1:3]
        self.step_keys = {}
        for ii in self._my_keys:
            self.step_keys[ii] = '--'.join(
                ["%s"%self.inputs.parameters["block_id"], ii]
            )
        
        self = _loop(
            self,
            self.step_keys,
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


class ConcurrentLearning(Steps):
    def __init__(
            self,
            name : str,
            block_op : Steps,
            step_config : dict = normalize_step_dict({}),
            upload_python_package : str = None,
    ):
        self.loop = ConcurrentLearningLoop(
            name+'-loop',
            block_op,
            step_config = step_config,
            upload_python_package = upload_python_package,
        )
        
        self._input_parameters={
            "type_map" : InputParameter(),
            "numb_models": InputParameter(type=int),
            "template_script" : InputParameter(),
            "train_config" : InputParameter(),
            "lmp_config" : InputParameter(),
            "fp_config" : InputParameter(),
        }
        self._input_artifacts={
            "exploration_scheduler" : InputArtifact(),
            "init_models" : InputArtifact(optional=True),
            "init_data" : InputArtifact(),
            "iter_data" : InputArtifact(),
            "fp_inputs" : InputArtifact(),
        }
        self._output_parameters={
        }
        self._output_artifacts={
            "exploration_scheduler": OutputArtifact(),
            "models": OutputArtifact(),
            "iter_data" : OutputArtifact(),
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

        self._init_keys = ['scheduler', 'id']
        self.loop_key = 'loop'
        self.step_keys = {}
        for ii in self._init_keys:
            self.step_keys[ii] = '--'.join(['init', ii])

        self = _dpgen(
            self,
            self.step_keys,
            name, 
            self.loop,
            self.loop_key,
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


def _loop (
        steps, 
        step_keys,
        name : str,
        block_op : OP,
        step_config : dict = normalize_step_dict({}),
        upload_python_package : str = None,
):    
    step_config = deepcopy(step_config)
    step_template_config = step_config.pop('template_config')
    step_executor = init_executor(step_config.pop('executor'))

    block_step = Step(
        name = name + '-block',
        template = block_op,
        parameters={
            "block_id" : steps.inputs.parameters["block_id"],
            "type_map" : steps.inputs.parameters["type_map"],
            "numb_models" : steps.inputs.parameters["numb_models"],
            "template_script" : steps.inputs.parameters["template_script"],
            "train_config" : steps.inputs.parameters["train_config"],
            "lmp_config" : steps.inputs.parameters["lmp_config"],
            "conf_selector" : steps.inputs.parameters["conf_selector"],
            "fp_config" : steps.inputs.parameters["fp_config"],
        },
        artifacts={
            "lmp_task_grp" : steps.inputs.artifacts["lmp_task_grp"],
            "fp_inputs" : steps.inputs.artifacts["fp_inputs"],            
            "init_models": steps.inputs.artifacts["init_models"],
            "init_data": steps.inputs.artifacts["init_data"],
            "iter_data": steps.inputs.artifacts["iter_data"],
        },
        key = step_keys['block'],
    )
    steps.add(block_step)

    scheduler_step = Step(
        name = name + '-scheduler',
        template=PythonOPTemplate(
            SchedulerWrapper,
            python_packages = upload_python_package,
            **step_template_config,
        ),
        parameters={
            "exploration_report": block_step.outputs.parameters['exploration_report'],
        },
        artifacts={
            "exploration_scheduler": steps.inputs.artifacts['exploration_scheduler'],
            "trajs" : block_step.outputs.artifacts['trajs'],
        },
        key = step_keys['scheduler'],
        executor = step_executor,
        **step_config,
    )
    steps.add(scheduler_step)

    id_step = Step(
        name = name + '-make-block-id',
        template=PythonOPTemplate(
            MakeBlockId,
            python_packages = upload_python_package,
            **step_template_config,
        ),
        parameters={
        },
        artifacts={
            "exploration_scheduler": scheduler_step.outputs.artifacts['exploration_scheduler'],
        },
        key = step_keys['id'],
        executor = step_executor,
        **step_config,
    )
    steps.add(id_step)

    next_step = Step(
        name = name+'-next',
        template = steps,
        parameters={
            "block_id" : id_step.outputs.parameters['block_id'],
            "type_map" : steps.inputs.parameters["type_map"],
            "numb_models" : steps.inputs.parameters["numb_models"],
            "template_script" : steps.inputs.parameters["template_script"],
            "train_config" : steps.inputs.parameters["train_config"],
            "lmp_config" : steps.inputs.parameters["lmp_config"],
            "conf_selector" : scheduler_step.outputs.parameters["conf_selector"],
            "fp_config" : steps.inputs.parameters["fp_config"],
        },
        artifacts={
            "exploration_scheduler" : scheduler_step.outputs.artifacts["exploration_scheduler"],
            "lmp_task_grp" : scheduler_step.outputs.artifacts["lmp_task_grp"],
            "fp_inputs" : steps.inputs.artifacts["fp_inputs"],
            "init_models" : block_step.outputs.artifacts['models'],
            "init_data" : steps.inputs.artifacts['init_data'],
            "iter_data" : block_step.outputs.artifacts['iter_data'],
        },
        when = "%s == false" % (scheduler_step.outputs.parameters['converged']),
    )
    steps.add(next_step)    

    steps.outputs.artifacts['exploration_scheduler'].from_expression = \
        if_expression(
            _if = (scheduler_step.outputs.parameters['converged'] == True),
            _then = scheduler_step.outputs.artifacts['exploration_scheduler'],
            _else = next_step.outputs.artifacts['exploration_scheduler'],
        )
    steps.outputs.artifacts['models'].from_expression = \
        if_expression(
            _if = (scheduler_step.outputs.parameters['converged'] == True),
            _then = block_step.outputs.artifacts['models'],
            _else = next_step.outputs.artifacts['models'],
        )
    steps.outputs.artifacts['iter_data'].from_expression = \
        if_expression(
            _if = (scheduler_step.outputs.parameters['converged'] == True),
            _then = block_step.outputs.artifacts['iter_data'],
            _else = next_step.outputs.artifacts['iter_data'],
        )

    return steps


def _dpgen(
        steps, 
        step_keys,
        name,
        loop_op,
        loop_key,
        step_config : dict = normalize_step_dict({}),
        upload_python_package : str = None
):    
    step_config = deepcopy(step_config)
    step_template_config = step_config.pop('template_config')
    step_executor = init_executor(step_config.pop('executor'))

    scheduler_step = Step(
        name = name + '-scheduler',
        template=PythonOPTemplate(
            SchedulerWrapper,
            python_packages = upload_python_package,
            **step_template_config,
        ),
        parameters={
            "exploration_report": None,
        },
        artifacts={
            "exploration_scheduler": steps.inputs.artifacts['exploration_scheduler'],
            "trajs" : None,
        },
        key = step_keys['scheduler'],
        executor = step_executor,
        **step_config,
    )
    steps.add(scheduler_step)

    id_step = Step(
        name = name + '-make-block-id',
        template=PythonOPTemplate(
            MakeBlockId,
            python_packages = upload_python_package,
            **step_template_config,
        ),
        parameters={
        },
        artifacts={
            "exploration_scheduler": scheduler_step.outputs.artifacts['exploration_scheduler'],
        },
        key = step_keys['id'],
        executor = step_executor,
        **step_config,
    )
    steps.add(id_step)

    loop_step = Step(
        name = name + '-loop',
        template = loop_op,
        parameters = {
            "block_id" : id_step.outputs.parameters['block_id'],
            "type_map" : steps.inputs.parameters['type_map'],
            "numb_models" : steps.inputs.parameters['numb_models'],
            "template_script" : steps.inputs.parameters['template_script'],
            "train_config" : steps.inputs.parameters['train_config'],
            "conf_selector" : scheduler_step.outputs.parameters['conf_selector'],
            "lmp_config" : steps.inputs.parameters['lmp_config'],
            "fp_config" : steps.inputs.parameters['fp_config'],
        },
        artifacts={
            "exploration_scheduler" : scheduler_step.outputs.artifacts['exploration_scheduler'],
            "lmp_task_grp" : scheduler_step.outputs.artifacts['lmp_task_grp'],
            "fp_inputs" : steps.inputs.artifacts['fp_inputs'],
            "init_models": steps.inputs.artifacts["init_models"],
            "init_data": steps.inputs.artifacts["init_data"],
            "iter_data": steps.inputs.artifacts["iter_data"],
        },
        key = '--'.join(["%s"%id_step.outputs.parameters['block_id'], loop_key]),
    )
    steps.add(loop_step)

    steps.outputs.artifacts["exploration_scheduler"]._from = \
        loop_step.outputs.artifacts["exploration_scheduler"]
    steps.outputs.artifacts["models"]._from = \
        loop_step.outputs.artifacts["models"]
    steps.outputs.artifacts["iter_data"]._from = \
        loop_step.outputs.artifacts["iter_data"]
    
    return steps