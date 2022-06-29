import json
import os
import sys
import numpy as np
from typing import List, Union, Tuple, Dict, Optional
from rid.utils import list_to_string
from rid.common.protein import get_dihedral_from_resid
from rid.common.plumed.plumed_constant import (
    dihedral_name,
    dihedral_def_from_atoms,
    deepfe_def,
    print_def,
    restraint_def,
    restraint_prefix
)

angle_id = {
    "phi": 0,
    "psi": 1
}

def make_restraint(
        name: str,
        arg: str,
        kappa: Union[str, int, float],
        at: Union[str, int, float]
    ):
    return restraint_def.format(
                name = name,
                arg = arg,
                kappa = kappa,
                at = at
            )

def make_restraint_list(
        cv_list: List[str],
        kappa: List[Union[int, float, str]],
        at: List[Union[int, float, str]]
    ) -> Tuple[List, List]:
    res_names = []
    res_list = []
    assert len(cv_list) == len(kappa), "Make sure `kappa` and `cv_names` have the same length."
    assert len(cv_list) == len(at), "Make sure `at` and `cv_names` have the same length."
    for idx, cv_print in enumerate(cv_list):
        res_name = "res-" + cv_print
        res_names.append(res_name)
        res_list.append(make_restraint(res_name, cv_print, kappa[idx], at[idx]))
    return res_list, res_names


def make_deepfe_bias(
        cv_list: List[str],
        trust_lvl_1: float = 1.0,
        trust_lvl_2: float = 2.0,
        model_list: List[str] = ["graph.pb"]
    ) -> str:
    if len(model_list) == 0:
        return ""
    cv_string = list_to_string(cv_list, ",")
    model_string = list_to_string(model_list, ",")
    return deepfe_def.format(
        trust_lvl_1 = trust_lvl_1,
        trust_lvl_2 = trust_lvl_2,
        model = model_string,
        arg = cv_string
    )


def make_print(
        name_list,
        stride,
        file_name
    ) -> str:
    return print_def.format(
        stride = stride,
        arg = list_to_string(name_list, ","),
        file = file_name
    )


def make_wholemolecules(atom_index):
    arg_list = list_to_string(atom_index, ",")
    return ("WHOLEMOLECULES" + " " +
            "ENTITY0=" + arg_list +
            "\n")


def user_plumed_def(cv_file, pstride, pfile):
    ret = ""
    cv_names = []
    print_content = None
    with open(cv_file, 'r') as fp:
        for line in fp.readlines():
            if ("PRINT" in line) and ("#" not in line):
                print_content = line + "\n"
                break
            ret += line + "\n"
            if (":" in line) and ("#" not in line):
                cv_names.append(("{}".format(line.split(":")[0])).strip())
            elif ("LABEL" in line) and ("#" not in line):
                cv_names.append(("{}".format(line.split("LABEL=")[1])).strip())
    if ret == "" or cv_names == "":
        raise RuntimeError("Invalid customed plumed files.")
    if print_content is not None:
        assert len(print_content.split(",")) == len(cv_names), "There are {} CVs defined in the plumed file, while {} CVs are printed.".format(len(cv_names), len(print_content.split(",")) )
        print_content_list = print_content.split()
        print_content_list[-1] = "FILE={}".format(pfile)
        print_content_list[1] = "STRIDE={}".format(str(pstride))
        print_content = " ".join(print_content_list)
    return ret, cv_names, print_content


def make_torsion(
        name: str,
        atom_list: List[Union[int, str]]
    ) -> str:
    assert len(atom_list) == 4, f"Make sure dihedral angle defined by 4 atoms, not {len(atom_list)}."
    return dihedral_def_from_atoms.format(
        name = name,
        a1 = atom_list[0], a2 = atom_list[1],
        a3 = atom_list[2], a4 = atom_list[3],
    )


def make_torsion_name(resid: int, angid: int):
    return dihedral_name.format(
        resid = resid,
        angid = angid
    )


def make_torsion_list(
        dihedral_info: Dict,
    ) -> Tuple[List, List]:
    torsion_list = []
    torsion_name_list = []
    for resid in dihedral_info.keys():
        for ang in dihedral_info[resid].keys():
            torsion_name = make_torsion_name(resid=resid, angid = angle_id[ang])
            torsion_name_list.append(torsion_name)
            torsion_list.append(make_torsion(
                name = torsion_name,
                atom_list=list(dihedral_info[resid][ang])
            ))
    return torsion_list, torsion_name_list


def make_torsion_list_from_file(
        file_path: str,
        selected_resid: List[int]
    ) -> Tuple[List, List]:
    return make_torsion_list(
        get_dihedral_from_resid( file_path, selected_resid
    ))


def make_restraint_plumed(
        conf: Optional[str] = None,
        cv_file: Optional[str] = None,
        selected_id: Optional[List[int]] = None,
        kappa: Union[int, float, List[Union[int, float]]] = 0.5,
        at: Union[int, float, List[Union[int, float]]] = 1.0,
        stride: int = 100,
        output: str = "plm.out",
        mode: str = "torsion"    
    ):
    content_list = []
    if mode == "torsion":
        cv_content_list, cv_name_list = \
            make_torsion_list_from_file(conf, selected_id)
        content_list += cv_content_list
    elif mode == "custom":
        ret, cv_name_list, _ = user_plumed_def(cv_file, stride, output)
        content_list.append(ret)
    else:
        raise RuntimeError("Unknown mode for making plumed files.")

    if not isinstance(kappa, List):
        kappa = [kappa for _ in range(len(cv_name_list))]
    if not isinstance(at, List):
        at = [at for _ in range(len(cv_name_list))]
    res_list, _ = make_restraint_list(
        cv_name_list, kappa, at
    )
    content_list += res_list
    content_list.append( make_print(cv_name_list, stride, output) )
    return list_to_string(content_list, split_sign="\n")


def make_deepfe_plumed(
        conf: Optional[str] = None,
        cv_file: Optional[str] = None,
        selected_resid: Optional[List[int]] = None,
        trust_lvl_1: float = 1.0,
        trust_lvl_2: float = 2.0,
        model_list: List[str] = ["graph.pb"],
        stride: int = 100,
        output: str = "plm.out",
        mode: str = "torsion"
    ):
    content_list = []
    if mode == "torsion":
        cv_content_list, cv_name_list = \
            make_torsion_list_from_file(conf, selected_resid)
        content_list += cv_content_list
    elif mode == "custom":
        ret, cv_name_list, _ = user_plumed_def(cv_file, stride, output)
        content_list.append(ret)
    else:
        raise RuntimeError("Unknown mode for making plumed files.")
    deepfe_string = make_deepfe_bias(cv_name_list, trust_lvl_1, trust_lvl_2, model_list)
    content_list.append(deepfe_string)
    content_list.append(make_print(cv_name_list, stride, output))
    return list_to_string(content_list, split_sign="\n")
