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
    mdp_content = "\n".join(content_list)
    return mdp_content


def modify_temperature(
        temperature: Union[str, float, int]
    ) -> Dict:
    return {
        "gen-temp": temperature,
        "ref-t": f"{temperature} {temperature}"
    }


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
        nsteps: Union[str, int], 
        output_freq: Union[str, int], 
        temperature: Union[str, int, float] = 300, 
        define: Optional[Union[str, List]] = None, 
        dt: float =0.002, 
        output_mode: str ="both"
    ):
    update_dict = {}
    update_dict.update({"nsteps": nsteps})
    update_dict.update({"dt": dt})
    if define is not None:
        update_dict.update(modify_define(define))
    update_dict.update(modify_temperature(temperature))
    update_dict.update(modify_output(output_freq, output_mode=output_mode))
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
    assert "temperature" in md_parameters_dict.keys()



if __name__ == "__main__":
    new = make_md_mdp_string(2000, 20)
    print(new)

