'''
Define mocked OP class from initial OP class
'''
import os
import shutil
from pathlib import Path
import json
from typing import Dict
import numpy as np
from dflow.python import (
    OP,
    OPIO,
    upload_packages
)
from dflow import (
    upload_artifact,
)

from rid.op.prep_data import CollectData, MergeData
from rid.utils import set_directory,save_txt
from rid.op.prep_exploration import PrepExplore
from rid.op.prep_label import PrepLabel
from rid.op.prep_label import CheckLabelInputs
from rid.op.run_exploration import RunExplore
from rid.op.run_label import RunLabel
from rid.op.label_stats import LabelStats
from rid.op.prep_select import PrepSelect
from rid.op.run_select import RunSelect
from rid.op.run_train import TrainModel
from rid.constants import (
    tf_model_name,
    init_conf_gmx_name,
    data_new,
    data_raw,
    gmx_conf_name,
    gmx_conf_out,
    plumed_output_name,
    gmx_mdrun_log,
    gmx_xtc_name,
    cv_force_out,
    mf_fig,
    cluster_fig,
    cluster_selection_data_name,
    cluster_selection_index_name,
    sel_gro_name,
    model_devi_name,
    sel_ndx_name,
    cv_init_label,
    model_devi_precision,
    train_log,
    train_fig,
    mf_std_fig
)

upload_packages.append(__file__)

def clear_dir(data=None):
    if data:
        if os.path.exists(data):
            shutil.rmtree(data)

def clear_files(files=None):
    if files:
        if isinstance(files,list):
            for file in files:
                if os.path.exists(file):
                    os.remove(file)
        elif isinstance(files,Path):
            if os.path.exists(files):
                    os.remove(files)

def make_mocked_init_models(numb_models):
    tmp_models = []
    for ii in range(numb_models):
        ff = Path(tf_model_name.format(tag=ii))
        ff.write_text(f'This is init model {ii}')
        tmp_models.append(ff)
    return tmp_models

def make_mocked_init_data(data_name,data_num):
    tmp_init_data = Path(data_name)
    tmp_init_data.mkdir(exist_ok=True, parents=True)
    tmp_datalist = []
    for ii in range(data_num):
        file_name = data_name+"_%s"%ii
        (tmp_init_data/file_name).write_text('this is %s %s'%(data_name,ii))
        tmp_datalist.append(tmp_init_data/file_name)
    
    return tmp_init_data,tmp_datalist

def make_mocked_init_confs(numb_confs):
    tmp_confs = []
    for ii in range(numb_confs):
        ff = Path(init_conf_gmx_name.format(idx=ii))
        ff.write_text(f'This is init conf {ii}')
        tmp_confs.append(ff)
    return tmp_confs

def make_mocked_file(file_name):
    ff = Path(file_name)
    ff.write_text(f'This is mocked %s'%file_name)
    return ff

def make_mocked_json(filename,fconts):
    ff = Path(filename)
    with open(ff,"w") as f:
        json.dump(fconts,f)
    return ff

class MockedCollectData(CollectData):
    @OP.exec_sign_check
    def execute(
            self,
            ip : OPIO,
    ) -> OPIO:
        with open(ip["cv_forces"][0], "r") as f:
            f1 = f.read()
        
        with open(data_new, "w") as f:
            f.write(f1)
    
        op = OPIO({
            "data_new" : Path(data_new)
            })
        return op

class MockedMergeData(MergeData):
    @OP.exec_sign_check
    def execute(
            self,
            ip : OPIO,
    ) -> OPIO:
        if ip["data_old"] is not None:
            with open(ip["data_old"], "r") as f:
                f1 = f.read()
        else:
            f1 = ""
        with open(ip["data_new"], "r") as f:
            f2 = f.read()
        
        with open(data_raw, "w") as f:
            if ip["data_old"] is not None:
                f.write(f1)
                f.write("\n")
            f.write(f2)
        
        op = OPIO({
            "data_raw" : Path(data_raw)
            })
        return op
    
class MockedPrepExplore(PrepExplore):
    @OP.exec_sign_check
    def execute(
            self,
            ip : OPIO,
    ) -> OPIO:
        task_path = Path(ip["task_name"])
        task_path.mkdir(exist_ok=True, parents=True)
        conf =  ip["conf"]
        with open(conf, "r") as f:
            fcons = f.read()
        with open(task_path.joinpath(gmx_conf_name), "w") as f:
            f.write(fcons)
        cv_dim = ip["cv_config"]["cv_dim"]
        op = OPIO({
                "task_path": task_path,
                "cv_dim": int(cv_dim)
            })
        print("finish prep!!!")
        return op

class MockedRunExplore(RunExplore):
    @OP.exec_sign_check
    def execute(
            self,
            ip : OPIO,
    ) -> OPIO:
        dp_conf_list = []
        assert(ip["task_path"].is_dir())
        print("start run explore !!")
        with set_directory(ip["task_path"]):
            with open(gmx_conf_name, "r") as f:
                fconts = f.read()
            print("write conf")
            with open(gmx_conf_out, "w") as f:
                f.write(fconts)
            print("write plumed")
            with open(plumed_output_name, "w") as f:
                f.write("plumed output")
            print("write log")
            with open(gmx_mdrun_log, "w") as f:
                f.write("gmx md log file")
            print("write traj")
            with open(gmx_xtc_name, "w") as f:
                f.write("gmx trajectory file")
        
        op = OPIO({
                "dp_selected_confs": dp_conf_list,
                "plm_out": ip["task_path"].joinpath(plumed_output_name),
                "md_log": ip["task_path"].joinpath(gmx_mdrun_log),
                "trajectory": ip["task_path"].joinpath(gmx_xtc_name),
                "conf_out": ip["task_path"].joinpath(gmx_conf_out)
            })
        return op

class MockedCheckLabelInputs(CheckLabelInputs):
    @OP.exec_sign_check
    def execute(
            self,
            ip : OPIO,
    ) -> OPIO:
        if ip["confs"] is None:
            if_continue = 0
            conf_tags = []
        else:
            if_continue = 1

            tags = {}
            for tag in ip["conf_tags"]:
                if isinstance(tag,Path):
                    with open(tag,"r") as f:
                        tags.update(json.load(f))
                else:
                    raise RuntimeError("Unkown Error.")
                
            conf_tags = []
            for conf in ip["confs"]:
                conf_tags.append(str(tags[conf.name]))
        
        op = OPIO(
            {
                "if_continue": if_continue,
                "conf_tags": conf_tags
            }
        )
        return op


class MockedPrepLabel(PrepLabel):
    @OP.exec_sign_check
    def execute(
            self,
            ip : OPIO,
    ) -> OPIO:
        assert(ip["at"].is_file())
        task_path = Path(ip["task_name"])
        task_path.mkdir(exist_ok=True, parents=True)
        conf =  ip["conf"]
        with open(conf, "r") as f:
            fcons = f.read()
        with open(task_path.joinpath(gmx_conf_name), "w") as f:
            f.write(fcons)
        op = OPIO({
                "task_path": task_path,
            })
        print("finish prep!!!")
        return op
    
class MockedRunLabel(RunLabel):
    @OP.exec_sign_check
    def execute(
            self,
            ip : OPIO,
    ) -> OPIO:
        assert(ip["task_path"].is_dir())
        print("start run label !!")
        with set_directory(ip["task_path"]):
            print("write plumed")
            with open(plumed_output_name, "w") as f:
                f.write("plumed output")
            print("write log")
            with open(gmx_mdrun_log, "w") as f:
                f.write("gmx md log file")
            with open(cv_force_out,"w") as f:
                f.write('1.000000e+00 2.000000e+00\n')
            with open("mf_info.out", "w") as f:
                f.write("mean force information file")
            with open(gmx_xtc_name, "w") as f:
                f.write("gmx trajectory")
        
        op = OPIO({
                "plm_out": ip["task_path"].joinpath(plumed_output_name),
                "trajectory": ip["task_path"].joinpath(gmx_xtc_name),
                "md_log": ip["task_path"].joinpath(gmx_mdrun_log),
                "cv_forces": ip["task_path"].joinpath(cv_force_out),
                "mf_fig": ip["task_path"].joinpath(mf_fig),
                "mf_info": ip["task_path"].joinpath("mf_info.out")
            })
        return op
    
class MockedLabelStats(LabelStats):
    @OP.exec_sign_check
    def execute(
            self,
            ip : OPIO,
    ) -> OPIO:
        cv_forces_list = [_ for _ in ip["cv_forces"] if _ is not None]
        task_path = Path("label_std")
        task_path.mkdir(exist_ok=True, parents=True)
        with set_directory(task_path):
            with open(mf_std_fig, "w") as f:
                f.write("mean force std fig")
        op = OPIO({
                "mf_std_fig":  task_path.joinpath(mf_std_fig),
                "cv_forces": list(cv_forces_list)  
            })
        return op
    
class MockedPrepSelect(PrepSelect):
    @OP.exec_sign_check
    def execute(
            self,
            ip : OPIO,
    ) -> OPIO:
        cluster_threshold = 2
        cls_sel_idx = np.array([0])
        numb_cluster = len(cls_sel_idx)
        selected_data = np.array([[1,2,3]])
        task_path = Path(ip["task_name"])
        task_path.mkdir(exist_ok=True, parents=True)
        with set_directory(task_path):
            with open(cluster_fig, "w") as f:
                f.write("this is the cluster fig")
            np.save(cluster_selection_index_name, cls_sel_idx)
            np.save(cluster_selection_data_name, selected_data)
        
        op_out = OPIO({
                "cluster_threshold": cluster_threshold,
                "numb_cluster": numb_cluster,
                "cluster_fig": task_path.joinpath(cluster_fig),
                "cluster_selection_index": task_path.joinpath(cluster_selection_index_name),
                "cluster_selection_data": task_path.joinpath(cluster_selection_data_name)
            })
        return op_out

class MockedRunSelect(RunSelect):
    @OP.exec_sign_check
    def execute(
            self,
            ip : OPIO,
    ) -> OPIO:
        task_path = Path(ip["task_name"])
        task_path.mkdir(exist_ok=True, parents=True)
        
        cls_sel_idx = np.load(ip["cluster_selection_index"])
        cls_sel_data = np.load(ip["cluster_selection_data"])
        
        sel_idx = cls_sel_idx
        sel_data = cls_sel_data
        conf_list = []
        cv_init_list = []
        conf_tags = {}
        with set_directory(task_path):
            save_txt(sel_ndx_name, sel_idx, fmt="%d")
            save_txt("cls_"+model_devi_name, [], fmt=model_devi_precision)
            for ii, sel in enumerate(sel_idx):
                make_mocked_file(sel_gro_name.format(walker = int(ip['task_name']),idx=sel))
                conf_list.append(task_path.joinpath(sel_gro_name.format(walker = int(ip['task_name']),idx=sel)))
                conf_tags[sel_gro_name.format(walker = int(ip['task_name']),idx=sel)] = f"{ip['task_name']}_{sel}"
                save_txt(cv_init_label.format(walker = int(ip['task_name']),idx=sel), sel_data[ii])
                cv_init_list.append(task_path.joinpath(cv_init_label.format(walker = int(ip['task_name']),idx=sel)))
        with open("conf.json", "w") as f:
            json.dump(conf_tags,f)

        op_out = OPIO(
            {
               "selected_confs": conf_list,
               "selected_cv_init": cv_init_list,
               "model_devi": task_path.joinpath("cls_"+model_devi_name),
               "selected_indices": task_path.joinpath(sel_ndx_name),
               "selected_conf_tags": task_path.joinpath("conf.json")
            }
        )
        return op_out

class MockedTrain(TrainModel):
    @OP.exec_sign_check
    def execute(
            self,
            ip : OPIO,
    ) -> OPIO:
        task_path = Path(ip["model_tag"])
        task_path.mkdir(exist_ok=True, parents=True)
        assert(ip["data"].is_file())
        with set_directory(task_path):
            tf_model = tf_model_name.format(tag=ip["model_tag"])
            with open(tf_model,"w") as f:
                f.write("this is trained model")
            train_log_name = train_log.format(tag=ip["model_tag"])
            with open(train_log_name,"w") as f:
                f.write("this is log file")
            train_fig_name = train_fig.format(tag=ip["model_tag"])
            with open(train_fig_name,"w") as f:
                f.write("this is train fig")
        
        op_out = OPIO(
            {
                "model": task_path.joinpath(tf_model),
                "train_log": task_path.joinpath(train_log_name),
                "train_fig": task_path.joinpath(train_fig_name)
            }
        )
        return op_out
