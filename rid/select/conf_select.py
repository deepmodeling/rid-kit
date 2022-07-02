import os, sys
import numpy as np
import logging
from typing import List, Union, Sequence
from rid.select.model_devi import make_std
from rid.constants import sel_gro_name
import mdtraj as md


logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


class ConfSelector:
    def __init__(
            self,
            threshold: float,
            model_list: List
        ):
        self.threshold = threshold
        self.model_list = model_list
    
    def select(
            self,
            data: Union[np.ndarray, List],
        ):
        stds = make_std(data, models=self.model_list)
        # selected_idx, selected_data = select_from_devi(data, stds, self.threshold)
        # return selected_idx, selected_data
    

def slice_xtc(
        xtc: str,
        top: str,
        selected_idx: Sequence
    ):
    traj = md.load_xtc(xtc, top=top)
    for sel in selected_idx:
        frame = traj[sel]
        frame.save_gro(sel_gro_name.format(idx=sel))


def select_from_devi(model_devi, threshold):
    selected_idx = []
    for idx, _std in enumerate(model_devi):
        if _std > threshold:
            selected_idx.append(idx)
    selected_idx = np.array(selected_idx, dtype=int)
    message = "max std: %f, min std: %f, avg std %f \n" % (
        np.max(model_devi), np.min(model_devi), np.average(model_devi))
    message += "number of angles than %f is %d" % (threshold, len(selected_idx))
    logger.info(message)
    return selected_idx


def adjust_trust_lvl(lvl1, lvl2, origin_lvl1, origin_lvl2, numb_cluster, threshold):
    # adaptive trust level
    if numb_cluster < threshold:
        enhc_trust_lvl_1 = lvl1 * 1.5
        enhc_trust_lvl_2 = lvl2 * 1.5
        # enhc_trust_lvl_2 = enhc_trust_lvl_1 + 1 
        # TODO: would it be OK to use 1?
    else:
        enhc_trust_lvl_1 = origin_lvl1
        enhc_trust_lvl_2 = origin_lvl2
    if enhc_trust_lvl_1 > origin_lvl1 * 8:
        enhc_trust_lvl_1 = origin_lvl1
        enhc_trust_lvl_2 = origin_lvl2
    return enhc_trust_lvl_1, enhc_trust_lvl_2

# def adjust_trust_lvl_mode(numb_cluster, threshold):
#     # adaptive trust level
#     mode = None
#     if numb_cluster < threshold:
#         return "increase"
#     else:
#         enhc_trust_lvl_1 = origin_lvl1
#         enhc_trust_lvl_2 = origin_lvl2
#     if enhc_trust_lvl_1 > origin_lvl1 * 8:
#         enhc_trust_lvl_1 = origin_lvl1
#         enhc_trust_lvl_2 = origin_lvl2
#     return enhc_trust_lvl_1, enhc_trust_lvl_2
