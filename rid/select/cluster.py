import os, sys
from typing import Union, List, Optional
import logging
import numpy as np
import sklearn.cluster as skcluster
from matplotlib import pyplot as plt
from rid.constants import cluster_fig
from pathlib import Path


logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


class Cluster:
    def __init__(
        self,
        cvs: Union[np.ndarray, List],
        threshold: float,
        task_name: str,
        angular_mask: Optional[Union[np.ndarray, List]] = None, 
        weights: Optional[Union[np.ndarray, List]] = None,
        max_search_step: int = 500,
        max_selection: int = 1000
    ):
        if angular_mask is None:
            angular_mask = np.zeros(shape=(cvs[1],))
        if weights is None:
            weights = np.ones(shape=(cvs[1],))
        self.angular_mask = angular_mask
        self.weights = weights
        self.task_name = task_name
        self.max_search_step = max_search_step
        self.threshold = threshold
        self.cvs = cvs
        self.enlarge_coeff = 1.05
        self.reduce_coeff = 0.95
        self.cls_sel = None
        self.max_selection = max_selection
    
    def make_threshold(self, numb_cluster_lower, numb_cluster_upper):
        current_iter = 0
        logger.info(f"set numb_cluster_upper to {numb_cluster_upper}")
        logger.info(f"set numb_cluster_lower to {numb_cluster_lower}")
        assert numb_cluster_lower < numb_cluster_upper, f"expect numb_cluster_upper > numb_cluster_lower, "
        "got {numb_cluster_upper} < {numb_cluster_lower}"
        while current_iter < self.max_search_step:
            logger.info(f"making threshold attempt {current_iter}")
            cls_sel = sel_from_cluster(
                self.cvs, self.threshold, Path(self.task_name),angular_mask=self.angular_mask, 
                weights=self.weights, max_selection=self.max_selection)
            test_numb_cluster = len(set(cls_sel))
            if test_numb_cluster < numb_cluster_lower:
                self.threshold = self.threshold * self.reduce_coeff
            elif test_numb_cluster > numb_cluster_upper:
                self.threshold = self.threshold * self.enlarge_coeff
            else:
                break
            logger.info(f"set threshold to {self.threshold}, get {test_numb_cluster} clusters.")
            current_iter += 1
        self.cls_sel = cls_sel
        return self.threshold
    
    def get_cluster_selection(self):
        if self.cls_sel is None:
            self.cls_sel = sel_from_cluster(
                self.cvs, self.threshold, Path(self.task_name),angular_mask=self.angular_mask, 
                weights=self.weights, max_selection=self.max_selection)
        return self.cls_sel
        

def cv_dist(cv1, cv2, angular_mask, weights):
    diff = cv1 - cv2
    angular_mask = np.array(angular_mask)
    angular_boolean = (angular_mask == 1)
    angular_diff = diff[angular_boolean]
    angular_diff[angular_diff < -np.pi] += 2 * np.pi
    angular_diff[angular_diff >  np.pi] -= 2 * np.pi
    diff[angular_boolean] = angular_diff
    return np.linalg.norm(diff * weights)


def mk_dist(cv, angular_mask, weights):
    nframe = cv.shape[0]
    dist = np.zeros([nframe, nframe])
    for ii in range(nframe):
        for jj in range(ii+1, nframe):
            dist[ii][jj] = cv_dist(cv[ii], cv[jj], angular_mask, weights)
            dist[jj][ii] = dist[ii][jj]
    return dist


def mk_cluster(dist, distance_threshold):
    logger.info("clustering ...")
    cluster = skcluster.AgglomerativeClustering(n_clusters=None,
                                          linkage='average',
                                          affinity='precomputed',
                                          distance_threshold=distance_threshold)
    cluster.fit(dist)
    return cluster.labels_


def chooseClusterCenter(dist:np.ndarray, conf_ids:list):
    id_min = conf_ids[0]
    min_loss = None
    # n_tot_frame = dist.shape[0]
    for id in conf_ids:
        rmsd = dist[id]
        loss = 0
        for i, conf in enumerate(conf_ids):
            rms = rmsd[conf]
            loss += rms * rms
        loss = np.sqrt(loss)
        if min_loss == None:
            min_loss = loss
        else:
            if loss < min_loss:
                id_min = id
                min_loss = loss
    return [id_min]


def sel_from_cluster(cvs, threshold, task_path, angular_mask=None, weights=None, max_selection=1000):
    if len(cvs) <= 1:
        return cvs
    weights = np.array(weights)
    dist = mk_dist(cvs, angular_mask, weights) 
    labels = mk_cluster(dist, threshold)
    # plot clustering distributions
    xlist = [i for i in range(len(labels))]
    plt.figure(figsize=(10, 8), dpi=100)
    plt.xlabel("trajectory frames")
    plt.ylabel("cluster index")
    plt.title("cluster distributions along trajectories")
    plt.scatter(xlist, labels, s = 5)
    plt.savefig(task_path.joinpath(cluster_fig))
    # make cluster map
    _cls_map = []
    for _ in range(len(set(labels))):
        _cls_map.append([])
    for idx, label in enumerate(labels):
        _cls_map[label].append(idx)
    cls_map = []
    for clust in _cls_map:
        cls_map.append((clust, len(clust)))
    cls_map = sorted(cls_map, key=lambda x: x[1], reverse=True)
    # randomly select from clusters
    cls_sel = []
    np.random.seed(seed=None)
    for cluster, _ in cls_map:
        # _ret = np.random.choice(cluster, 1)
        _ret = chooseClusterCenter(dist, cluster)
        cls_sel.append(_ret[0])
    if len(cls_sel) > max_selection:
        cls_sel = cls_sel[:max_selection]
        logger.info("selection number is beyond max selection, adjust to the max number.")
    return np.array(cls_sel, dtype=int)
