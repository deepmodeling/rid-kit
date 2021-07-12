# RiD Tutorial

In this tutorial, we will go over how RiD runs iteratively. Every procedure of rid will be introduced, including how it works, what input it needs, what it will output.

## Preparation

To begin with, we need to prepare the following things:
> 1. Initial conformation files in the PDB format.
> 2. A RiD configuration file, offering parameters used in the whole procedure.
> 3. A `cv.json` file containing the Collective Variables (CVs) selected.
> 4. A machine configuration file, setting the localhost or dispatching systems.

### Global Setting in RiD

RiD needs to update its NN iteratively, thus, we need to set the number of iteration in `rid.json` by adjusting the key `numb_iter`.

When first sampling (e.g. module enhcMD), several independent and individual trajectories will go from the same or different initial conformations. These initial conformation files in PDB format should be put into a directory (e.g. `rid-kit/examples/mol` where we offer 3 same ala2 conformation files as instances). To specify the number of individual trajectories which we call as a 'walker', we can modify the key `numb_walkers` in `rid.json: "00.bias"`. NOTE: the number of initial conformation files should not be less than `numb_walkers`. If these walkers start from the same origin conformation, just make enough copies of the initial conformation file. Then, set the number of neural networks by `"numb_model"`, and set number of nodes of each layer by `"neurons"` like `"neurons": [ 200, 200, 200, 200 ]`. Other parameters are also available in `rid.json`.

### Machine Setting

We will introduce other parameters in `rid.json` in the following sections. But now, we have to turn to `machine.json` or `local.json` to prepare for our starting. Actually, these two files are the same essentially. But we still give these two demos. Each of them looks like:
```JSON
    "JobName": {
        "machine":{
            "batch_type": "Slurm",
            "context_type": "LazyLocalContext",
            "local_root": "./",
            "remote_root": "./"
        },
        "resources":{
            "queue_name": "GPU_2080Ti",
            "number_node": 1,
            "cpu_per_node": 8,
            "gpu_per_node": 1,
            "group_size": 1,
            "if_cuda_multi_devices": false
        }
    },
```
For every job we run, we can assign a host for it. `"batch_type"` is the type of your dispatching system (e.g. `Slurm` or `PBS`), or if you run jobs on localhost, set `"batch_type": "Shell"`.

`resources` represents the resources you apply from the dispatching system. However, if you run them on localhost, parameters `queue_name`, `number_node`, `cpu_per_node`, `gpu_per_node` are meaningless, because, `dpdispatcher` will occupy full resources for these jobs by default. In this case, just set `queue_name` to any `str`, and set the following three parameters to any `int` number.

RiD will distribute jobs by the package `dpdispatching`. In each task (e.g. `enhcMD`), several jobs (e.g. `gmx mdrun`) will be dispatched. `dpdispatching` will group these jobs and submit these reformed job groups to compute nodes or localhost. `group_size` means how many jobs contained per group. If you set `group_size` to a very large number, it will be the same as a sequential run. 

### Selecting Collective Variables (CVs)

In this version, the users can choose the dihedral angles as CVs.  In the CV file(`cv.json`), the users can write the indexes of the selected residues, the two dihedral angles ($\psi$ and $\phi$) will be chosen as the CVs. 

Let's begin with a simple example, ala2, which has a sequence (1ACE, 2ALA, 3NME). The `cv.json` file can be set as:
```JSON
{
    "_comment":	   " dihedral angles: phi, psi ",
    "dih_angles" : [ [ {"name" : ["C"],  "resid_shift" : -1},
		       {"name" : ["N"],  "resid_shift" : 0},
		       {"name" : ["CA"], "resid_shift" : 0},
		       {"name" : ["C"],  "resid_shift" : 0} ], 
		     [ {"name" : ["N"],  "resid_shift" : 0}, 
		       {"name" : ["CA"], "resid_shift" : 0},
		       {"name" : ["C"],  "resid_shift" : 0},
		       {"name" : ["N"],  "resid_shift" : 1} ]
		   ],
    "selected_index":  [0, 1, 2],
    "alpha_idx_fmt":	"%03d",
    "angle_idx_fmt":	"%02d"
}
```
`"dih_angles"` is our definition of dihedral angles($\phi$, $\psi$) by default. Users can write the list of `"selected_index"` as their wish. Rid-kit will remove the non-existed dihedral angles of the terminal residues automatically. In this example, `"selected_index":  [0, 1, 2]` means we select dihedral angles of the 1st, 2nd and 3rd residues as our CVs. However, the terminal residues (or caps) have only either $\phi$ or $\psi$, or none of them (e.g. 1ACE and 3NME have no dihedral angles, 2ALA has $\phi$ and $\psi$), so even if we have selected the indexes of 1ACE and 3NME, the total dimension of CVs is **2**, which comes from the two dihedral angles of 2ALA.  

> ***Note***: The indexes in `cv.json` start from **0**, while the indexes of residues in `.gro` file start from **1**.

## Get Started
We need to import all the main blocks we will need. If you have installed successfully, the following imports should never fail. 
`from rid import enhcMD, resMD, train, gen_rid`

Also, we need some auxiliary functions, which will read checkpoint, record tasks, and help build the format.
`from rid.lib.utils import record_iter, get_checkpoint, make_iter_name`

Some other third-party packages should be available.
```python
import json, os, glob
import argparse
```
Then, set a valid path for your outputs. In the examples we offered, this variable is named to `out_dir` which will be the basepath for all files.

`record_file` will record the checkpoint(iteration indexes and task indexes) of tasks after finishing one. If you want to rerun the process and make sure that a record file exists in the work path, the program will restart from the next one of the end of the record(just use the same command to restart). If a task was restarted, but a working folder (which this task should generate) has already existed, this existed folder will be backed off as `folder_name.bk00`. That is, you can restart the process at any individual task node by modifying the recording file.

However, if there is NOT a record file in the working path, the whole process will restrat at the very beginning. The old one will become a backup folder as `old_folder.bk000`.


## Main Procedure

> ***Note***: 
> 
> We would still recommend you go through these articles before starting this section.
> 
>  [1]  Zhang, L., Wang, H., E, W.. Reinforced dynamics for enhanced sampling in large atomic and molecular systems[J]. The Journal of chemical physics, 2018, 148(12): 124113.
>  
>  [2]  Wang, D., Zhang, L., Wang, H., E, W.. Efficient sampling of high-dimensional free energy landscapes using adaptive reinforced dynamics[J]. arXiv preprint arXiv:2104.01620, 2021.
>

After having files above well-prepared, we can now perform rid.
```python
from rid.lib.utils import record_iter, get_checkpoint, make_iter_name
```
These are useful tools for us to get checkpoints and formulate our folder names.

Main modules should also be imported.
```python
from rid import enhcMD, resMD, train, gen_rid
```
`enhcMD` is the module of enhanced sampling; `resMD` is the module of calculating the free energy; `train` is the module of training the neural networks. `gen_rid` is to generate the working folder.

As a examples, we set our variables as:
```python
out_dir    # the output path.
mol_dir    # the folder containing initial conformation files.
rid_json   # the rid configuration file with the JSON format.
cv_file    # the cv file with the JSON format.
```

then, get the checkpoint:
```python
record_file = os.path.join(out_dir, record_name)  # record_name is the name of recording file, "record.txt" by default.
checkpoint = get_checkpoint(record_file)  # check point in list format. It will be set a negative number if we first strat our jobs.

if sum(checkpoint) < 0:  # if true, it means we start from the beginning.
    print("prepare gen_rid")
    gen_rid(out_dir, mol_dir, rid_json)  # check the conformation files and generate the working folder.
```

Since we don't have any neural networks in the first interation, so the `init_model` is empty.
```python
prev_model = init_model
```

Now, we can design a circulation:
```python
max_tasks = 10
number_tasks = 8  # by default, you can adjust this number according to your demand.
for iter_idx in range(iter_numb):
    # we need to get the neural networks we have trained after the 1st iteration. This code help get the .pb file generated after training.
    if iter_idx > 0 :
        prev_model = glob.glob(out_dir + "/" + make_iter_name(iter_idx-1) + "/02.train/*pb")
    
    # We just 
    for tag in range(number_tasks):
        if iter_idx * max_tasks + tag <= checkpoint[0] * max_tasks + checkpoint[1]:
            continue  # if true, this task has finished.
        elif tag == 0:
            ...
        elif tag == 1:
            ...
```

### enhcMD

Then, we can arrange our tasks in order.
```python
#  Following codes above

elif tag == 0:  # prepare
    print("prepare gen_enhc")
    enhcMD.make_enhc(iter_idx, rid_json, prev_model, mol_dir, cv_file ,base_dir=out_dir)
    record_iter(record_file, iter_idx, tag)
elif tag == 1:  # md run
    print("run enhanced MD")
    enhcMD.run_enhc(iter_idx, rid_json, machine_json, base_dir=out_dir)
    record_iter(record_file, iter_idx, tag)
elif tag == 2:  # post process
    print("prepare post enhc")
    enhcMD.post_enhc(iter_idx, rid_json, machine_json, base_dir=out_dir)
    record_iter(record_file, iter_idx, tag)
```
This part (3 steps) help us sample(sampling with brute force in the 1st iteration and enhanced sampling after the 1st sampling) and split the trajectory into separated .gro files. The sub working folder is named as `00.enhcMD` under the main working folder and there will be folder named from `000` representing each walker under the sub-folder. `"bias_nsteps"` and `"bias_frame_freq"` help us set the total steps and the output frequency in the bias MD.

Several files will be generated. Let's see files under `(output)/00.enhcMD/000/` for example. You can find a folder named `confs` easily which contains all the separated .gro file for each frame in the trajectory. `.log` file record all the output of Gromacs. `plm.out` is the output of Plumed and we will use this file as our training data.

If you run this part in the second iteration, you could see `.pb` files (soft-link) under the path like `graph.000.pb`. That means we have used the neural networks that we have tarined in the previous iterations, so the process would be enhanced. Remember that we have a variable defined as `prev_model = glob.glob (out_dir + "/" + make_iter_name(iter_idx-1) + "/02.train/*pb")`, and this could help us find the graph files we need.

We follow the steps of adpative RiD, so the trust level are adjusted adaptively. `"bias_trust_lvl_1"` and `"bias_trust_lvl_2"` are the initial guesses of trust level. After every sampling, we will cluster the conformations we have selected. If the cluster number is less than a thrashold(`"num_of_cluster_threshold"` in `rid.json`), the trust level will be adjusted. We will use agglomerative clustering algorithm to cluster, which will be introduced in the next section. Fianlly, the trust level that has been adjusted will be recorded in the file `trust_lvl1.dat`.


### resMD

Next step, we will select conformations and calculate the free energy of the conformations we select. Here comes `resMD`.

```python
#  Following codes above

 elif tag == 3:  # preparation
    print("prepare gen_res")
    resMD.make_res(iter_index=iter_idx, json_file=rid_json, cv_file=cv_file, mol_path=mol_dir, base_dir=out_dir)
    record_iter(record_file, iter_idx, tag)
elif tag == 4:  # calculation of free energy
    print("prepare run res")
    resMD.run_res (iter_index=iter_idx, json_file=rid_json, machine_json=machine_json, base_dir=out_dir)
    record_iter(record_file, iter_idx, tag)
elif tag == 5:  # post process
    print("prepare post res")
    resMD.post_res(iter_index=iter_idx, json_file=rid_json, machine_json=machine_json, cv_file=cv_file, base_dir=out_dir)
    record_iter(record_file, iter_idx, tag)
```

In the step of preparation, we will select conformations sampled in `enhcMD` that fall into the areas with high in NN. In the first iteration, all conformations may be selected since anywhere on the free erergy surface has high uncertainty. However, in the following iterations, we will only select conformations falling into areas with high-uncertainty by calculating the standard deviation of the neural networks (usually 4 neural networks). The threshold in this selection process is set as `"sel_threshold"` in `rid.json`. Any structure (or sampled CV) having higher standard deviation than `"sel_threshold"` will be selected. Still under the path `(work_path)/00.enhcMD/00x/`,  `sel.log` records the process of selection. `"sel.out"` records the indexes of frames that we select. 

After that, we will cluster the conformations that we have selected. You could set the initial interval of the cluster number in `rid.json` for the first clustering. `"init_numb_cluster_upper"` and `"init_numb_cluster_lower"` will give us a interval like `[init_numb_cluster_lower, init_numb_cluster_upper]`, and rid will adjust the threshold automatically to let the number of clusters in the first iteration fall into the inteval. This threshold will be recorded as `cluster_threshold.dat` in the path `(work_path)/00.enhcMD/00x/` and used in the following iterations.  Note that, there is another value in `rid.json` also named as `"cluster_threshold"`, but this value is the initial guess and will be adjusted according to the cluster number and the interval you set. `cls.sel.out` is final results that containing the conformtaions you have selected randomly from every cluster, and RiD will select one conformation from each cluster. `resMD.make_res` will finish these procedures.

Then we will use restrained MD to calculate the free energy of the conformations according to the indexes in `sel.out`. You will see many sub-folders under `work_path/01.resMD/` has the format `xxx.xxxxxx`. The number in front of the dot indicates the walker index, and the number after the dot represents the index of the frame. `resMD.run_res` helps us finish them. All parameters(steps, time interval, force constant, etc.) about MD could be set in `rid.json`.

Change path into one of them, you will see all the files about resMD. Actually, we need mean force as our final training data. `force.out` in each folder records this value. `centers.out` records the center of the CVs. `kappa.out` records the force constant. It is clear that diviation from the center times the force constant is the mean force. `resMD.post_res` will finish these procedures.

After processing, One could see `.../01.resMD/data.raw` easily, which is the training data we will use in the next section.

### Train

In this section, we will handle our neural networks.

```python
#  Following codes above

elif tag == 6:
    print("preparation of training")
    train.make_train(iter_index=iter_idx, json_file=rid_json, base_dir=out_dir)
    record_iter(record_file, iter_idx, tag)
elif tag == 7:
    print("run training")
    train.run_train(iter_index=iter_idx, json_file=rid_json, machine_json=machine_json, cv_file=cv_file, base_dir=out_dir)
    record_iter(record_file, iter_idx, tag)
```

Remember that, at the beginning of the tutorial, we have set `"numb_model"` and `"neurons"`. You see some sub folders representing every model. There is also a sub folder named as `data` containing the data, in which `data.new.raw` represents the data selected in the current iteration, `data.old.raw` contains the data selected in previous iterations and `data.raw` are the mixed file of new data and old data. Note that, there is a key `"res_iter"` in `rid.json`. If the current index of iteration is less than this value, we will simply mix all of the data and use the original parameters you have set; but if the current index of iteration is greater than this value, we will only use part of the old data, and this ratio can be set by `"res_olddata_ratio"` which means `(old data used) : (old data not used) = (res_olddata_ratio) : 1`. 

`train.make_train` helps us mix data and process checkpoint files. `train.run_train` will train the networks according to the
parameters. Finally, you can see several model files named as `graph.xxx.pb`. These files will be copied to next iteration to add biased potential. `train.log` will record the error during the training process.


Now, all steps are done. 







