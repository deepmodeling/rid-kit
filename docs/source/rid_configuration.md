# Setting Parameters of `RiD-kit`

`RiD-kit` uses a JSON-format file (typically `rid.json`) to configure simulations. Here we explain these parameters one by one.

## Overall Settings

* **`name`** `(str)` is the name of this task. 

* **`numb_walkers`** `(int)` is the number of parallel walkers to explore. `Exploration` Step can be achieved by `numb_walkers` parallel trajectories simultaneously. Data of these parallel walkers are collected together into to train neural networks.

* **`numb_iters`** `(int)` is the maximum number of iterations of RiD workflow. As it is very convenient to continue and rerun the RiD workflow, it does not really matter to set a accurate value. A recommended value is greater than 5, typically around 10 for the first attemption.

* **`trust_lvl_1`** `(int)` and **`trust_lvl_2`** `(int)` (or $e_0$ and $e_1$ in published papers) are two thresholds to control the biased forces and select data. In biased simulation, the bias forces are tuned by model deviations:
  $$
  F(r) = -\nabla_{r_i} U(r) + \sigma( \epsilon ( s( r))) \nabla_{r_i} A(r) \\
  \sigma(\epsilon)=
            \begin{cases}
                    1, & \epsilon<\epsilon_0 \\
                    \frac{1}{2}+\frac{1}{2}\cos{(\pi \frac{\epsilon-\epsilon_0}{\epsilon_1-\epsilon_0})}, & \epsilon_0 <\epsilon < \epsilon_1 \\
                    0, &\epsilon > \epsilon_1
            \end{cases}
  $$
  In data selection, data will be collected if their model deviations are greater than `trust_lvl_1`.

  In adaptive RiD version, these two values refer to the initial trust levels and will be adjusted according to the number of clusters during simulations.

* **`init_models`** `(List[str])` are the initial guesses of neural networks. Usually we know nothing about the systems and `[]` is set to it.

### Example

```JSON
"name": "test",
"numb_walkers": 2,
"numb_iters": 20,
"trust_lvl_1": 2,
"trust_lvl_2": 3,
"init_models": [],
```

## CV setting

This section configures collective variables (CVs). `RiD-kit` provides three modes to configure CV: `"torsion"`, `"distance"` and `"custom"`. 

### Torsion mode

In torsion mode, `RiD-kit` uses torsion (dihedral angles) of proteins as collective variables. Set `"mode": "torsion"`. `selected_resid`, `angular_mask` and `weights` must not be `none` or empty if you use torsion mode.

* **`selected_resid`** `(List[int])` residue ids (starting form 1) of selected residues. Two dihedral angles of each selected residue, $\phi$ and $\psi$, are used. Note that the first residue of a chain (N terminal) doesn't have $\phi$ and the last residue of a chain (C terminal) doesn't have $\psi$.

* **`angular_mask`** `(List[int])` the mask of augular (periodic) CVs, 1 for periodic and 0 for non-periodic. In torsion mode, all CVs are periodic, so a list filled by 1 with length equal to number of CVs should be set.

* **`weights`** `(List[int])` weights of CVs to scale their values. Used in clustering to calculate the Euclidean distances between CVs. This can prevent from CV discrimination if some CV's range is much larger than another one.

### Distance mode
In distance mode, `RiD-kit` uses distance between atoms of systems as collective variables. Set `"mode": "distance"`. 
`selected_atomid`, `angular_mask` and `weights` must not be `none` or empty if you use distance mode.

* **`selected_atomid`** `(List[List[int]])` ids (starting form 1) of selected pair of atoms. 

* **`angular_mask`** `(List[int])` the mask of augular (periodic) CVs, 1 for periodic and 0 for non-periodic. In distance mode, all CVs are nonperiodic, so a list filled by 0 with length equal to number of CVs should be set.

* **`weights`** `(List[int])` weights of CVs to scale their values. Used in clustering to calculate the Euclidean distances between CVs. This can prevent from CV discrimination if some CV's range is much larger than another one.

### Custom mode

`RiD-kit` also supports user-defined collective variables.
Set `"mode": "custom"` to use customed CV of your own design. `cv_file`, `angular_mask` and `weights` must not be `none` or empty if you use custom mode. The CV file is in the PLUMED2 format, you should add your own CV to the `PRINT` line in the CV file.

Note that if you use only one CV file, the CV file name should not end with `".pdb"`. If you use multiple files to define your own CV, the file to define your CV should not end with `".pdb"`, and other files should end with `".pdb"`.

* **`"cv_file"`** `(List[str])` List of Paths to CV files. The files define collective variables in PLUMED2 format. Technically, CV that PLUMED2 supports can be suporrted by `RiD-kit`. 
  
* **`"angular_mask"`** `(List[int])` the mask of augular (periodic) CVs, 1 for periodic and 0 for non-periodic. In custom mode, figure out the periodic CVs and set 1 at the corresponding location in list.
  
* **`"weights"`** `(List[int])` the same as above.
 
### Example
```JSON
"CV": {
    "mode": "torsion",
    "selected_resid": [ 1, 2 ],
    "angular_mask": [ 1, 1 ],
    "weights": [ 1, 1 ],
    "cv_file":[""]
}
```
```JSON
"CV": {
        "mode": "distance",
        "selected_atomid": [[2,5],[5,7]],
        "angular_mask": [0,0],
        "weights": [1,1],
        "cv_file":[""]
}
```
```JSON
"CV": {
        "mode": "custom",
        "selected_atomid":[[161,165],[124,156],[124,161],[156,165]],
        "angular_mask": [0,0,0,0],
        "weights": [1,1,1,1],
        "cv_file": ["colvar", "plmpath.pdb"]
}
```

## ExploreMDConfig
This section configures the parameters in `Exploration` step. Currently, rid-kit supports two types of sampler: `"gmx"` and `"lmp"`, stand for Gromacs and Lammps respectively.

### gmx type

Set `"type": "gmx"` to use Gromacs as sampler. `nstep`, `temperature`, `ref-t`,`output_freq`, `dt`, `output_mode` must not be `none` or empty if you use gmx type. You can also set other parameters in `mdp` files for your usage (Acutally it is recommended since the default settings in rid-kit may be be the best suit for your system).

* **`nsteps`** `(int)` Number of steps of MD simulation in exploration step.

* **`type`** `(str)` Type of sampler in `Exploration` step. Currently rid-kit support `gmx` (Gromacs) and `lmp` (Lammps).

* **`temperature`** `(int)` Temperature of MD simulations. Please make sure this value is the same as the temperature in `Label` step unless you are meant to keep them different.

* **`ref-t`** `(str)` Temperature for each group of the system, default is `300 300` for `water non-water` group. 

* **`output_freq`** `(int)`  Frame output frequence of MD simulations. A recommended value is `nsteps/1000` to make sure at least 1000 frames generated during exploration.

* **`dt`** `(int)` Time interval of MD simulations in `ps`. 0.002 is recommended for normal simulations. One may use a larger interval, e.g. 0.004, when heavy hydrogen modes in Gromacs.

* **`output_mode`** `(str):` Optional modes: `"both", "single", "double", "none"`. 
  * `"both"`: Generate both full presicion format `.trr` and compressed presision format `.xtc` trajectories during MD simulations.
  * `"single"` Only generate compressed presision format `.xtc` trajectories during MD simulations.
  * `"double"` Only generate full presicion format `.trr` trajectories during MD simulations.
  *  `"none"` Don't generate trajectory files. (Used for tasks that only need PLUMED2 ouput.)
  
* **`ntmpi`** `(int)` Number of thread-MPI ranks to start (0 is guess). See detail in [Gromacs manual](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html).
* **`nt`** `(int)` Total number of threads to start (0 is guess). See detail in [Gromacs manual](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html).
* **`max_warning`** Max warnings in `gmx grompp` steps. See detail in [Gromacs manual](https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html).

### lmp type
Set `"type": "lmp"` to use Lammps as sampler. `output_freq`, `dt` and `inputfile` must not be `none` or empty if you use lmp type. 

* **`inputfile`** `(str)` Name of lammps input file to define the simulation paramers.

### Example

```JSON
"ExploreMDConfig": {
        "nsteps": 50000,
        "type": "gmx",
        "temperature": 300,
        "output_freq": 50,
        "ref-t": "300 300",
        "verlet-buffer-tolerance":"-1",
        "rlist": 1,
        "rvdw": 0.9,
        "rvdw-switch": 0,
        "rcoulomb": 0.9,
        "rcoulomb-switch": 0,
        "epsilon-r":1,
        "epsilon-rf":80,
        "dt": 0.002,
        "fourierspacing": "0.12",
        "output_mode": "single",
        "ntmpi": 1,
        "nt": 8,
        "max_warning": 2
}
```
```JSON
"ExploreMDConfig": {
        "type":"lmp",
        "inputfile":"input_explore.lammps",
        "dt": 0.002,
        "output_freq": 50,
        "ntmpi": 1,
        "nt": 8,
        "max_warning": 2
}
```

## SelectorConfig
This section configures the parameters in `Selection` step. In `Selection` step, all CV values are clustered. Then data owning high model deviation are selected, collected and sent to `Label` step. 

* **`"cluster_threshold"`** `(int)` Initial guess of cluster threshold. Note: the real cluster threshold is generated from this guess.

* **`numb_cluster_lower`** `(int)` and **`numb_cluster_upper`** `(int)` These two values form an closed interval `[numb_cluster_lower, numb_cluster_upper]` to make a proper cluster threshold. From the initial guess of cluster, threshold will be adjusted to let the number of clusters fall into the interval. This process only happens in the first iteration. The threshold will be fixed in the following iterations where the `trust level` will be adjusted  in adaptive version of RiD. See [published paper](https://doi.org/10.1038/s43588-021-00173-1) for detail.

* **`"max_selection"`** `(int)` The max selection number of clusters during `Selection` step. If number of clusters is greater than this threshold, the first `max_selection`th clusters will be selected.

* **`numb_cluster_threshold`** `(int)` If number of clusters of MD trajectories in exploration step at current interation is less than this value, the trust level will be adjusted. See [published paper](https://doi.org/10.1038/s43588-021-00173-1) for detail. A recommended value is half of `numb_cluster_lower`.

* **`slice_mode`** `(str)` Optional values: `"gmx"` and `"mdtraj"`. `RiD-kit` extracts selected frame from MD trajectorie. `gmx` mode uses, Gromacs `gmx trjconv` to slice trajectories, `mdtraj` mode uses `mdtraj` python interface to slice trajectories. We highly recommed using `gmx` mode due to known bugs ([#Issue1514 ](https://github.com/mdtraj/mdtraj/issues/1514)) from `mdtraj` of changing `.gro` topology names. 


### Example

```JSON
"SelectorConfig": {
    "cluster_threshold": 1,
    "numb_cluster_lower": 16,
    "numb_cluster_upper": 26,
    "max_selection": 30,
    "numb_cluster_threshold": 8,
    "slice_mode": "gmx"
},
```

## LabelMDConfig
This section configures the parameters in `Label` step. Currently, rid-kit supports two methods of labeling: `"restrained"` and `"constrained"`, stand for `restrained MD` and `constrained MD` respectively. Most settings are quite similar to those in `Exploration` Step. The simulation time is usually different between `Exploration` and `Label` steps, while the simulation time for `Exploration` has more freedom, the simulation time for `Label` step has to be chosen with care (longer enough to ensure convergence and shorter enough to avoid wasting). Our experience is that `100ps` is enough for `torsion mode` in `restrained MD method`, `1ns` is enough for `distance mode` in `constrained MD method`.  

### restrained method
Set `"method": "restrained"` to use restrained MD as mean force calculator. The only different parameters with `Exploration` step is `kappas` and `std_threshold`.

* **`kappas`** `(List[int])` A list of force constants ($\kappa$) of harmonic restraints. The length of the list is equal to the number of CVs.

* **`std_threshold`** `(float)`(default 5.0, the unit is consistent with the mean force) A number represents the mean force standard deviation threshold, beyond which the mean force is neglected and will not be used in the dataset for training free energy model. You should test labeling MD for your own system to determine an appropriate number for this threshold.

### constrained method
Set `"method": "constrained"` to use constrained MD as mean force calculator. Currently rid-kit only supports distance CV to use this method, also only `gmx` type is supported to perform constrained MD simulation. The other parameters is the same with `Exploration` step.

Note that if you want to use constrained MD as the mean force calculator, apart from setting `method` to be `constrained` in the `label_config`, you should add `[ constraints ]` line corresponding to the `[ moleculartype ]` in your input `topology` file yourself, since gromacs specifies constraints information for each `[ moleculartype ]`.

Also note that, since gromacs only supports constrained MD for `distance CV`, the constrained MD simulation in rid-kit only supports `distance CV` at this moment.


### Example

```JSON
    "LabelMDConfig": {
        "nsteps": 50000,
        "temperature":300,
        "method": "restrained",
        "type": "gmx",
        "output_freq": 100,
        "ref-t": "300 300",
        "rlist": 1,
        "verlet-buffer-tolerance":"-1",
        "rvdw": 0.9,
        "rvdw-switch": 0,
        "rcoulomb": 0.9,
        "rcoulomb-switch": 0,
        "epsilon-r":1,
        "epsilon-rf":80,
        "dt": 0.002,
        "fourierspacing": "0.12",
        "output_mode": "single",
        "ntmpi": 1,
        "nt": 8,
        "max_warning": 2,
        "kappas": [ 500, 500 ],
        "std_threshold": 2.0
    }
```
```JSON
"LabelMDConfig": {
        "nsteps": 50000,
        "temperature":300,
        "method": "constrained",
        "type": "gmx",
        "output_freq": 100,
        "ref-t": "300 300",
        "rlist": 1,
        "verlet-buffer-tolerance":"-1",
        "rvdw": 0.9,
        "rvdw-switch": 0,
        "rcoulomb": 0.9,
        "rcoulomb-switch": 0,
        "epsilon-r":1,
        "epsilon-rf":80,
        "dt": 0.001,
        "fourierspacing": "0.12",
        "output_mode": "both",
        "ntmpi": 1,
        "nt": 8,
        "max_warning": 2,
        "std_threshold": 10.0
}
```

## Train
This section configures the parameters in `Train` step. `RiD-kit` is based on `Tensorflow`.

* **`numb_models`** `(int)` Number of models that are trained in `Train` step. `RiD-kit` uses model deviations (or standrad deviation of output of these models) to evaluate the quality of free energy surface, so `numb_models` mush be **greater** than 1.

* **`neurons`** `(List[int])` The number of neurons of each layer. `RiD-kit` uses MLP as the basic neural network structure. Number of elements in list means the number of hidden layers and each element defines number of nodes in each layer. For example, `[ 50, 50, 50, 50 ]` means there are 4 hidden layers and each hidden layers has 50 neurons.

* **`resnet`** `(bool)` Wether to use residual connection between layers. If `true`, the number of nodes of layers must be equal.

* **`epoches`** `(int)` Numebr of epoches.

* **`init_lr`** `(float)` Initial learning rate. It will decay exponentially during training.

* **`decay_steps`** `(int)` Decay steps of learning rate. See [tensorflow api docs](https://www.tensorflow.org/api_docs/python/tf/compat/v1/train/exponential_decay) for detail.

* **`decay_rate`** `(float)` Decay rate of learning rate. See [tensorflow api docs](https://www.tensorflow.org/api_docs/python/tf/compat/v1/train/exponential_decay) for detail.

* **`drop_out_rate`** `(float)` Dropout rate of dropout layers.

* **`numb_threads`** `(int)` Threads of training.

### Example

```JSON
"Train": {
    "numb_models": 4,
    "neurons": [ 50, 50, 50, 50 ],
    "resnet": true,
    "batch_size": 32,
    "epoches": 2000,
    "init_lr": 0.0008,
    "decay_steps": 120,
    "decay_rate": 0.96,
    "drop_out_rate": 0.1,
    "numb_threads": 8,
    "use_mix": false,
    "restart": false
}
```

## Full Examples of `rid.json`

You can find full examples of `rid.json` within `"rid-kit/rid/template"`.