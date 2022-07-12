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


def select_from_devi(model_devi, threshold):
    logger.info("select conformations by model deviations.")
    selected_idx = []
    for idx, _std in enumerate(model_devi):
        if _std > threshold:
            selected_idx.append(idx)
    selected_idx = np.array(selected_idx, dtype=int)
    logger.info("max std: %f, min std: %f, avg std %f \n" % (
        np.max(model_devi), np.min(model_devi), np.average(model_devi)))
    logger.info("number of angles than %f is %d" % (threshold, len(selected_idx)))
    return selected_idx
