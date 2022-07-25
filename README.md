
# Rid-kit
## 0. **Content**
<!-- vscode-markdown-toc -->
* 1. [Introduction](#Introduction)
* 2. [Installation](#Installation)
	* 2.1. [Environment installation](#Environmentinstallation)
        * 2.1.1. [Install python and tensorflow](#Installpythonandtensorflow)
        * 2.1.2. [Install tensorflow's C++ interface](#Installpythonandtensorflow)
        * 2.1.3. [Install plumed2.8.0](#Installplumed2.8.0)
        * 2.1.4. [Install gromacs 2021.4](#Installgromacs2021.4)
    * 2.2. [Install dpdispatcher](#Installdpdispatcher)
	* 2.3. [Install rid package](#Installridpackage)
* 3. [Quick Start](#QuickStart)
	* 3.1. [CV selection](#CVselection)
	* 3.2. [Dispatching](#Dispatching)
		* 3.2.1. [Batch Job](#BatchJob)
		* 3.2.2. [Local Job](#LocalJob)
* 4. [Main procedure of RiD](#MainprocedureofRiD)
		* 4.1. [a. Biased MD](#a.BiasedMD)
		* 4.2. [b. Restrained MD](#b.RestrainedMD)
		* 4.3. [c. Neural network training](#c.Neuralnetworktraining)
* 5. [RiD settings](#RiDsettings)
	* 5.1. [rid.json](#rid.json)


##  1. <a name='Introduction'></a>**Introduction**

Rid-kit is a python package for enhanced sampling via the RiD (Reinforced Dynamics) method.

##  2. <a name='Installation'></a>**Installation**

###  2.1. <a name='Environmentinstallation'></a>**Environment installation**

> ***Note***:
> 
> If you want to set environments by hand, please follow settings 1.1.1 ~ 1.1.4. 
You can also install environment of rid through the dockerfile, just run "docker build --net=host --tag rid-kit ." and then run "docker run -d -it --gpus all --net=host -v /mnt/workdir/:/mnt/workdir/ --shm-size=400g --name rid-container rid-kit" you can get the required enviroment.
#### 2.1.1 <a name='Installpythonandtensorflow'></a>**Install python and tensorflow**

#### 2.1.2 <a name='Install tensorflows C++ interface'></a>**Install tensorflow's C++ interface**
##### 2.1.2.1 <a name='Using Bazel'></a>**Using Bazel**
It is relatively simple and one can follow the steps <https://docs.deepmodeling.com/projects/deepmd/en/master/install/install-tf.2.3.html>

The problem with Bazel compile is that it requires accessing website aboard. One way out is to use the compiled file and set path variables `PATH` and `LD_LIBRARY_PATH` accordingly.

##### 2.1.2.2 <a name='Using Conda'></a>**Using Conda**
Use the command
```
conda install libtensorflow_cc=2.6.2=*cuda113* nccl -c conda-forge
```
to get the compiled libtensorflow_cc. Or one can create a conda enviroment 
```
conda create -n rid-drop python=3.9 libtensorflow_cc=2.6.2=*cuda110* nccl MDAnalysis mdtraj numpy
```
After installation set the ```tf_path``` as the ```CONDA_PREFIX```
```
conda activate rid-drop
export tf_path=${CONDA_PREFIX}
```
**You may encounter the error when compiling the plumed**
```
undefined reference to `std::__throw_bad_array_new_length()@glibcxx_3.4.29'
```
The reason for this error is that the tensorflow_cc library from conda requires GLIBCXX 3.4.29, but the highest version of GLICBCXX associated with the docker image is 3.4.28 (located in /usr/lib/x86_64-linux-gnu), so we need to upgrade the libstdc++6 library
```
DEBIAN_FRONTEND=noninteractive apt-get install -y software-properties-common
add-apt-repository -y ppa:ubuntu-toolchain-r/test
apt-get update
apt-get install -y libstdc++6
```

**You may encounter the error when compiling gromacs**
```
/opt/conda/lib/libcurl.so.4: no version information available (required by /usr/bin/cmake)
...
/opt/conda/lib/libtinfo.so.6: no version information available (required by /bin/bash)
```
The reason for this error is that there is version conflict of libcurl.so.4 and libtinfo.so.6 between the conda libraries and the libraries associated with ubuntu, we can relink the two libraries
```
rm ${CONDA_PREFIX}/lib/libtinfo.so.6
ln -s /usr/lib/x86_64-linux-gnu/libtinfo.so.6 ${CONDA_PREFIX}/lib/libtinfo.so.6
rm ${CONDA_PREFIX}/lib/libcurl.so.4
ln -s /usr/lib/x86_64-linux-gnu/libcurl.so.4 ${CONDA_PREFIX}/lib/libcurl.so.4
```






#### 2.1.3 <a name='Installplumed2.8.0'></a>**Install plumed2.8.0**
You need to copy compiled `DeePFE.cpp` to the plumed directory. This file locates at `rid-kit/install/DeePFE.cpp`
```bash
tar -xvzf plumed-2.8.0.tgz
cp DeePFE.cpp plumed-2.8.0/src/bias
cd plumed-2.8.0
./configure --prefix=$CONDA_PREFIX \
                   --enable-cxx=14 \
                   CXXFLAGS="-std=gnu++14 -I $CONDA_PREFIX/include/" \
                   LDFLAGS=" -L$CONDA_PREFIX/lib -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=$CONDA_PREFIX/lib/" \
make -j 6
make install
```
Set the ~/.bashrc
```bash
export LIBRARY_PATH=/opt/conda/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/conda/lib:$LD_LIBRARY_PATH
export PLUMED_KERNEL=/opt/conda/lib/libplumedKernel.@SOEXT@
```

#### 2.1.4 <a name='Installgromacs2021.4'></a>**Install gromacs 2021.4**

```bash
tar -xzvf gromacs-2021.4.tar.gz
cd gromacs-2021.4
plumed patch -p  # make sure you use the correct plumed version
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON \
                -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
                -DGMX_GPU=CUDA  \
                -DGMX_SIMD=avx_512
make -j 4
make install
```
Set the bashrc
```bash
source $CONDA_PREFIX/bin/GMXRC
```

### 2.2. <a name='Installdpdispatcher'></a>**Install dpdispatcher**
```bash
git clone https://github.com/deepmodeling/dpdispatcher.git
cd dpdispatcher
python setup.py install
```
dpdispatcher is a tool for job submission.

### 2.3. <a name='Installridpackage'></a>**Install rid package**
Now you have all dependence of RiD (Gromacs, Tensorflow, and a conda environment).
~~~bash
cd rit-kit
python setup.py install
~~~
Open python, try `import rid`.

Installation finishes successfully if you get no error.

##  3. <a name='QuickStart'></a>**Quick Start**
We offer a simple but complete example in `rid-kit/examples`. RiD process can run either on batch or locally.

if running locally, please try:
```bash
cd examples
python main.py jsons/rid.json -c jsons/cv.json -s jsons/local.json -i ./mol -o ./test_examples 
```
if running on batch, please try:
```bash
cd examples
python main.py jsons/rid.json -c jsons/cv.json -s jsons/machine.json -i ./mol -o ./test_examples 
```

To begin with, you should offer a rid parameters file(rid.json), a CV file(cv.json), a machine configuration file(machine.json), and a folder(mol/) containing initial conformation files in detail, and the number of conformation files should be equal to the number of walkers for parallel.

All these files are presented in `examples` folder where the users can adjust the parameter as their will.

The process will be recorded in `(work_path)/recoord.txt` in which its iteration indexes and task indexes will be written after finishing every task. If you want to rerun the process and make sure that a record file exists in the work path, the program will restart from the next one of the end of the record(just use the same command to restart). If a task was restarted, but a working folder (which this task should generate) has already existed, this existed folder will be backed off as `folder_name.bk00`. That is, you can restart the process at any individual task node by modifying the recording file.

However, if there is NOT a record file in the working path, the whole process will restrat at the very beginning. The old one will become a backup folder as `old_folder.bk000`.

###  3.1. <a name='CVselection'></a>**CV selection**
In this version, the users can choose the dihedral angles as CVs.  In the CV file(`cv.json`), the users can write the indexes of the selected residues, the two dihedral angles ($\psi$ and $\phi$) will be set as the CVs. 

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

Plumed will output all selected angles during every MD process, and the users can find them in `work_path/iter.0000xx/00.enhcMD/00x/plm.out`, file `angle.rad.out` in the same path is a copy but removing the frame indexes. Thus, in the previous example of ala2, the processed output `angle.rad.out` will look like:
```
-2.429137 2.552929
-2.503469 2.463779
...
-1.240340 2.390756
```
These data are nothing but the dihedral angles in every frame. The first column is $\phi$ angle, the second column is $\psi$ angle. 

We will add more features for users to select more different (and customed) CVs.

###  3.2. <a name='Dispatching'></a>**Dispatching**
Every task of RiD can be assigned to either compute nodes or local machines, which can be achieved in the machine configuration file(for instance, `machine.json or local.json`). These settings give the messages of users to `dpdispatcher` which can automatically distribute resources.

####  3.2.1. <a name='BatchJob'></a>**Batch Job**
If you want to submit jobs to a dispatching system, like Slurm or PBS, please follow settings like this:
```JSON
    "enhcMD": {
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
`"enhcMD"` represents the name of your jobs. Key `"machine"` decides where the jobs run, as for batch jobs, set `batch_type=your_dispatching_system`. Key `"resources"` means the resource you apply from the dispatching system. `"group_size"` is the size of each group of tasks which you have submitted, for example, you have 20 Tasks(a command, a Task, like `"gmx mdrun"`) and have `group_size = 10`, then the program will submit `20/10=2` Jobs to the dispatching system for parallel running, and each Job contains 10 Tasks(commands) going sequentially. At last, `"if_cuda_multi_devices"` can help assign tasks to one GPU.

####  3.2.2. <a name='LocalJob'></a>**Local Job**
If you want to run jobs locally, please follow settings like this:
```JSON
    "cmpf": {
        "machine":{
            "batch_type": "Shell",
            "context_type": "LazyLocalContext",
            "local_root": "./",
            "remote_root": "./"
        },
        "resources":{
            "queue_name": "queue_name",
            "number_node": 0,
            "cpu_per_node": 0,
            "gpu_per_node": 0,
            "group_size": 1000
        }
    },
```
Just set `batch_type` to `shell`. When you submit jobs to local machine, some keys(`"queue_name", "number_node", "cpu_per_node"`) become unnecessary. In these cases, the local machine will use full resources to carry out commands by default. If you have a quite large `group_size`, it will be the same as you run the command in a new shell.

By adjusting these settings, the users can assign any task anywhere.

##  4. <a name='MainprocedureofRiD'></a>**Main procedure of RiD**

RiD will run in iterations. Every iteration contains tasks below:

1. Biased MD;
2. Restrained MD;
3. Training neural network.

####  4.1. <a name='a.BiasedMD'></a>a. **Biased MD**

Just like Metadynamics, RiD will sample based on a bias potential given by NN models. An uncertainty indicator will direct the process of adding bias potential.

####  4.2. <a name='b.RestrainedMD'></a>b. **Restrained MD**

This procedure will calculate the mean force based on the sampling results, which can generate data set for training. 

####  4.3. <a name='c.Neuralnetworktraining'></a>c. **Neural network training**

A fully connected NN will be trained via sampling data. This network will generate a map from selected CV to free energy.

A more detailed description of RiD is published now, please see:

>  [1]  Zhang, L., Wang, H., E, W.. Reinforced dynamics for enhanced sampling in large atomic and molecular systems[J]. The Journal of chemical physics, 2018, 148(12): 124113.
>  
>  [2]  Wang, D., Zhang, L., Wang, H., E, W.. Efficient sampling of high-dimensional free energy landscapes using adaptive reinforced dynamics[J]. arXiv preprint arXiv:2104.01620, 2021.


##  5. <a name='RiDsettings'></a>**RiD settings**


Two necessary JSON files are required to get start a RiD procedure.

1. rid.json for configuration of simulation.
2. cv.json for specifying CV.

###  5.1. <a name='rid.json'></a>**rid.json**

**General setting**

| Parameters | Type | Description | Default/Example |
| :----: | :----: | :----: | :----: |
| gmx_prep | str | Gromacs preparation command | gmx grompp -maxwarn 1 |
| gmx_run | str | Gromacs md run command | gmx mdrun -ntmpi 1 |
| init_graph | list&str | initial graph files list | [] |
| numb_iter | int | number of iterations | 3 |
| numb_walkers | int | number of walkers | 2 |
| bf_traj_stride | int | brute force trajectory stride | 500 |

**Setting for biased MD**

| Parameters | Type | Description | Default/Example |
| :----: | :----: | :----: | :----: |
| bias_trust_lvl_1 | int | trust upper level | 2 |
| bias_trust_lvl_2 | int | trust lower level | 3 |
| bias_nsteps | int | total number of steps of biased MD | 20000 |
| bias_frame_freq | int | frame frequency for recording | 20 |
| sel_threshold | float/int | initial threshold for selection | 2 |
| cluster_threshold | float | * | 1.5 |
| num_of_cluster_threshold | int | minimum of cluster number | 8 |
| max_sel | int | maximum of selection of clusters | 30 |
| bias_dt | float | time interval of biased MD | 0.002 |
| bias_temperature | float/int | temperature for biased MD | 320 |

**Setting for restrained MD**

| Parameters | Type | Description | Default/Example |
| :----: | :----: | :----: | :----: |
| res_nsteps | int | total number of steps of restrained MD | 25000 |
| res_frame_freq | int | frame frequency for recording| 50 |
| res_dt | float | time interval of restrained MD | 0.002 |
| res_temperature | int | temperature for restrained MD | 320 |
| res_kappa | float/int | force constant for restraining CV | 500 |
| **res_traj_stride** | int | brute force trajectory stride | 500 |
| res_ang_stride | int | step stride of angle | 5 |
| res_prt_file | str | file name | plm.res.out |
| init_numb_cluster_upper | int | upper limit of cluster selection | 26 |
| init_numb_cluster_lower | int | lower limit of cluster selection, must be > 1 | 16 |
| conf_start | int | the index of the first conformation selected | 0 |
| conf_every | int | the stride of conformation selection | 1 |

**Setting for training and neural network**

| Parameters | Type | Description | Default/Example |
| :----: | :----: | :----: | :----: |
| numb_model | int | number of nn models | 4 |
| neurons | list&int | number of nodes for each layer | [256, 128, 64, 32] |
| resnet | bool | whether to use Resnet | True |
| batch_size | int | batch size | 128 |
| numb_epoches | int | total number of epochs for every training | 2000 |
| starter_lr | float | initial learning rate | 0.0008 |
| decay_steps | int | decay steps of lr | 120 |
| decay_rate | float | decay rate of lr | 0.96 |
| **res_iter** | int | after this iteration, old data will be reduced | 13 |
| res_numb_epoches | int | restrat setting | 2000 |
| res_starter_lr | float | restrat setting | 0.0008 |
| res_olddata_ratio | int/float | restrat setting | 7 |
| res_decay_steps | int | restrat setting | 120 |
| res_decay_rate | float | restrat setting | 0.96 |

