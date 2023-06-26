# Downloading files in `RiD-kit`

The download command is executed through
```bash
rid download WORKFLOW_ID -p STEP-KEY -f FILE_NAME -a ITERATION_START -e ITERATION_END -o OUTPUT_DIR
```

`ITERATION_START`(default 1) specifies from which iteration to start download, `ITERATION_END`(default 100) specifies to which iteration to end download. `OUTPUT_DIR`(default "./") is the output directories of the download.

The `STEP-KEY` in rid includes the following steps: `prep-exploration`, `run-exploration`, `prep-select`, `run-select`, `prep-label`, `run-label`, `label-stats`, `collect-data`, `merge-data`, `train`, `model-devi`.

The `FILE_NAME` in `prep-exploration` includes `task_path`. Which represents the input directory for Exploration MD simulation.

The `FILE_NAME` in `run-exploration` includes `trajectory`, `conf_out`,`plm_out`,`md_log`,`bias_fig`, `model_devi_fig`, `dp_model_devi_fig (optional)`, `dp_model_devi (optional)`, `dp_selected_indices (optional)`, `dp_selected_confs (optional)`, `projected_fig (optional)`.

The `trajectory` represents the trajectories of each walker produced by the Exploration MD, `conf_out` represents the last frame in the trajectories and will be used in the next iteration, `plm_out` represents the CV output from PLUMED, `md_log` is the log file from the molecular dynamics simulation, `bias_fig` is the figure showing the biasing potential used in the simulation, `model_devi_fig` is the figure showing the deviation of the free energy model used in the simulation, `dp_model_devi_fig (optional)` is the figure showing the deviation of the Deep Potential (DP) model used in the simulation, `dp_model_devi (optional)` contains numerical data on the deviation of the Deep Potential model used in the simulation, `dp_selected_indices (optional)` contains the indices of configurations with high dp model deviation selected for further analysis or training, `dp_selected_confs (optional)` contains the actual configurations high dp model deviation, `projected_fig (optional)` is the figure showing some projection of the simulation trajectory, currently only support for projection on the path CV. 

The `FILE_NAME` in `prep-select` includes `cluster_fig`, `cluster_selection_index`, `cluster_selection_data`.

The `cluster_fig`  is the figure showing the distribution of the clusters on the trajectories, `cluster_selection_index` contains the indices of the configurations from each cluster, `cluster_selection_data` contains the actual configurations selected from each cluster.

The `FILE_NAME` in `run-select` includes `selected_confs`, `selected_cv_init`, `model_devi`, `selected_indices`,`selected_conf_tags`.

The `selected_confs` contains the configurations that have been selected for labeling MD in the next step, `selected_cv_init` contains the initial values of the collective variables (CVs) for the selected configurations, `model_devi` contains the deviation of the free energy model on the selected configurations from the last `prep-select` step, `selected_indices` contains the indices of the selected configurations, `selected_conf_tags`  contains tags for the selected configurations. 

The `FILE_NAME` in `prep-label` includes `task_path`. Which represents the input directory for Exploration MD simulation.

The `FILE_NAME` in `run-label` includes `plm_out`, `cv_forces`, `mf_info`, `mf_fig`,`md_log`.

The `plm_out` represents the CV output from PLUMED, `cv_forces` contains the collective variables (CVs) and the mean forces obtained from the simulation, `mf_info` contains information about the mean force (MF) calculated during the simulation, including the values and the standard deviation, `mf_fig` is the figure showing the running average of the mean force as a function of time, `md_log`  is a log file from the MD simulation.

The `FILE_NAME` in `label-stats` includes `mf_std_fig` and `cv_forces`.

The `mf_std_fig` is a histogram showing the distribution of the standard deviation of the mean force, `cv_forces` is a filtered version of the same file in last `run-label` step, removing mean forces with high standard deviation.

The `FILE_NAME` in `collect-data` includes `data_new`. Which is the .npy file combining all the `cv_forces` files in the last step, represents the new training data.

The `FILE_NAME` in `merge-data` includes `data_raw`. Which is the .npy file containing all the training data in the previous iterations.

The `FILE_NAME` in `train` includes `model`, `train_log`, `train_fig`.

The `model` is the trained free energy model, `train_log` is the log of the training process. `train_fig` represents the training error during the training process.

The `FILE_NAME` in `model_devi` includes `model_devi`, `model_devi_fig`.

The `model_devi` is the model deviation on the trajectories of each walker using the newly trained free energy model, `model_devi_fig` is the figure of the model deviation.