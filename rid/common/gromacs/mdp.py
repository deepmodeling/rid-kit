from turtle import update
from rid.common.gromacs.gmx_constant import mdp_parameters
from typing import Optional, Dict, List, Union
    

def make_mdp_from_json(
        task: str, 
        inputs:Optional[Dict] = None
    ) -> Dict:
    assert task in mdp_parameters.keys()
    mdp_json = mdp_parameters[task]
    if inputs is not None:
        mdp_json.update(inputs)
    def make_mdp_line(key, value):
        return "{} \t= {}".format(key, value)
    content_list = []
    for key, value in mdp_json.items():
        content_list.append( make_mdp_line(key, value) )
    content_list.sort()
    mdp_content = "\n".join(content_list)
    return mdp_content

def modify_output(
        freq: Union[str, int], 
        output_mode: str = "both"
    ) -> Dict:
    if output_mode == "both":
        output_json = {
            "nstxout": freq,
            "nstvout": freq,
            "nstfout": freq,
            "nstenergy": freq,
            "nstxtcout": freq
        }
    elif output_mode == "single":
        output_json = {
            "nstxout": 0,
            "nstvout": 0,
            "nstfout": 0,
            "nstenergy": freq,
            "nstxtcout": freq
        }
    elif output_mode == "double":
        output_json = {
            "nstxout": freq,
            "nstvout": freq,
            "nstfout": freq,
            "nstenergy": freq,
            "nstxtcout": 0
        }
    elif output_mode == "none":
        output_json = {
            "nstxout": 0,
            "nstvout": 0,
            "nstfout": 0,
            "nstenergy": 0,
            "nstxtcout": 0
        }
    else:
        raise RuntimeError("Unknown output mode. Please specify one from 'single', 'double' or 'both'.")
    return output_json


def modify_define(
        define: Union[str, List]
    ) -> Dict:
    if type(define) == List:
        define_string = " ".join(define)
    else:
        define_string = define
    return {
        "define": define_string
    }

def make_md_mdp_string(
        gmx_config
    ):
    update_dict = {}
    for item in gmx_config:
        if item not in ["method","temperature","nt", "kappas","ntmpi", "max_warning", "output_freq", "output_mode","type","dp_model"]:
            update_dict[item] = gmx_config[item]
    update_dict.update(modify_output(freq = gmx_config["output_freq"], output_mode=gmx_config["output_mode"]))
            
    mdp_string = make_mdp_from_json(task="md", inputs=update_dict)
    return mdp_string


def make_md_mdp_from_config(
        md_parameters_dict
    ):
    check_basic_argument(md_parameters_dict)
    mdp_string = make_mdp_from_json(task="md", inputs=md_parameters_dict)
    return mdp_string


def check_basic_argument(md_parameters_dict: Dict):
    assert "nsteps" in md_parameters_dict.keys()
    assert "output_freq" in md_parameters_dict.keys()
