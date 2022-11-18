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

This section configures collective variables (CVs). `RiD-kit` provides two modes to configure CV: `"torsion"` and `"custom"`. 

### Torsion mode

In torsion mode, `RiD-kit` uses torsion (dihedral angles) of proteins as collective variables. Set `"mode": "torsion"`. `selected_resid`, `angular_mask` and `weights` must not be `none` or empty if you use torsion mode.

* **`selected_resid`** `(List[int])` residue ids (starting form 1) of selected residues. Two dihedral angles of each selected residue, $\phi$ and $\psi$, are used. Note that the first residue of a chain (N terminal) doesn't have $\phi$ and the last residue of a chain (C terminal) doesn't have $\psi$.

* **`angular_mask`** `(List[int])` the mask of augular (periodic) CVs, 1 for periodic and 0 for non-periodic. In torsion mode, all CVs are periodic, so a list filled by 1 with length equal to number of CVs should be set.

* **`weights`** `(List[int])` weights of CVs to scale their values. Used in clustering to calculate the Euclidean distances between CVs. This can prevent from CV discrimination if some CV's range is much larger than another one.

### Custom mode

`RiD-kit` also supports user-defined collective variables.
Set `"mode": "custom"`. `cv_file`, `angular_mask` and `weights` must not be `none` or empty if you use custom mode.

* **`"cv_file"`** `(str)` Path to a CV file. This file defines collective variables in PLUMED2 format. Technically, CV that PLUMED2 supports can be suporrted by `RiD-kit`. 
  
* **`"angular_mask"`** `(List[int])` the mask of augular (periodic) CVs, 1 for periodic and 0 for non-periodic. In custom mode, figure out the periodic CVs and set 1 at the corresponding location in list.
  
* **`"weights"`** `(List[int])` the same as above.
 
### Example
```JSON
"CV": {
    "mode": "torsion",
    "selected_resid": [ 1, 2 ],
    "angular_mask": [ 1, 1 ],
    "weights": [ 1, 1 ],
    "cv_file":""
},
```

## ExploreMDConfig

This section configures the parameters in `Exploration` step. MD eigine is Gromacs patched by PLUMED2. Parameter convention follows Gromacs and PLUMED2.


* **`nsteps`** `(int)` Number of steps of MD simulation in exploration step.

* **`output_freq`** `(int)`  Frame output frequence of MD simulations. A recommended value is `nsteps/1000` to make sure at least 1000 frames generated during exploration.

* **`temperature`** `(int)` Temperature of MD simulations. Please make sure this value is the same as the temperature in `Label` step unless you are meant to keep them different.

* **`dt`** `(int)` Time interval of MD simulations. 0.002 is recommended for normal simulations. One may use a larger interval, e.g. 0.004, when heavy hydrogen modes in Gromacs.

* **`output_mode`** `(str):` Optional modes: `"both", "single", "double", "none"`. 
  * `"both"`: Generate both full presicion format `.trr` and compressed presision format `.xtc` trajectories during MD simulations.
  * `"single"` Only generate compressed presision format `.xtc` trajectories during MD simulations.
  * `"double"` Only generate full presicion format `.trr` trajectories during MD simulations.
  *  `"none"` Don't generate trajectory files. (Used for tasks that only need PLUMED2 ouput.)
  
* **`ntmpi`** `(int)` Number of thread-MPI ranks to start (0 is guess). See detail in [Gromacs manual](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html).
* **`nt`** `(int)` Total number of threads to start (0 is guess). See detail in [Gromacs manual](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html).
* **`max_warning`** Max warnings in `gmx grompp` steps. See detail in [Gromacs manual](https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html).

### Example

```JSON
"ExploreMDConfig": {
    "nsteps": 25000,
    "output_freq": 25,
    "temperature": 300,
    "dt": 0.002,
    "output_mode": "both",
    "ntmpi": 1,
    "nt": 8,
    "max_warning": 0
},
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

This section configures the parameters in `Label` step. Most settings are quite similar to those in `Exploration` Step. In `Label`, `RiD-kit` performs restrained MD simulations, where hamonic restraints are exerted on CVs. These procedues need much shorter steps than `Exploration`.

* **`kappas`** `(List[int])` A list of force constants ($\kappa$) of harmonic restraints. The length of the list is equal to the number of CVs.

* **`nsteps`** `(int)` Number of steps of MD simulation of restrained simulations.

* **`output_freq`** `(int)`  Frame output frequence of MD simulations. A recommended value is `nsteps/1000` to make sure at least 1000 frames generated during exploration.

* **`temperature`** `(int)` Temperature of MD simulations. Please make sure this value is the same as the temperature in `Exploration` step unless you are meant to keep them different.

* **`dt`** `(int)` Time interval of MD simulations. 0.002 is recommended for normal simulations. Tune it down if you use a very large harmonic force, otherwise numerical explosion may occur.

* **`output_mode`** `(str):` Optional modes: `"both", "single", "double", "none"`. 
  * `"both"`: Generate both full presicion format `.trr` and compressed presision format `.xtc` trajectories during MD simulations.
  * `"single"` Only generate compressed presision format `.xtc` trajectories during MD simulations.
  * `"double"` Only generate full presicion format `.trr` trajectories during MD simulations.
  *  `"none"` Don't generate trajectory files. (Used for tasks that only need PLUMED2 ouput.)
  
* **`ntmpi`** `(int)` Number of thread-MPI ranks to start (0 is guess). See detail in [Gromacs manual](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html).
* **`nt`** `(int)` Total number of threads to start (0 is guess). See detail in [Gromacs manual](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html).
* **`max_warning`** Max warnings in `gmx grompp` steps. See detail in [Gromacs manual](https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html).


### Example

```JSON
"LabelMDConfig": {
    "kappas": [ 500, 500 ],
    "nsteps": 25000,
    "output_freq": 50,
    "temperature": 300,
    "dt": 0.002,
    "output_mode": "both",
    "ntmpi": 1,
    "nt": 8,
    "max_warning": 0
},
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

## A full Example of `rid.json`

You can find a full example of `rid.json` within `"rid-kit/rid/template"`. Or you can copy one from following:

```JSON
{
    "name": "test",
    "numb_walkers": 2,
    "numb_iters": 20,
    "trust_lvl_1": 2,
    "trust_lvl_2": 3,
    "init_models": [],
    
    "CV": {
        "mode": "torsion",
        "selected_resid": [ 1, 2 ],
        "angular_mask": [ 1, 1 ],
        "weights": [ 1, 1 ],
        "cv_file":""
    },

    "ExploreMDConfig": {
        "nsteps": 25000,
        "output_freq": 25,
        "temperature": 300,
        "dt": 0.002,
        "output_mode": "both",
        "ntmpi": 1,
        "nt": 8,
        "max_warning": 0
    },

    "SelectorConfig": {
        "numb_cluster_lower": 16,
        "numb_cluster_upper": 26,
        "cluster_threshold": 1,
        "max_selection": 30,
        "numb_cluster_threshold": 8,
        "slice_mode": "gmx"
    },

    "LabelMDConfig": {
        "nsteps": 25000,
        "output_freq": 50,
        "temperature": 300,
        "dt": 0.002,
        "output_mode": "both",
        "ntmpi": 1,
        "nt": 8,
        "max_warning": 0,
        "kappas": [ 500, 500 ]
    },

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
        "restart": false,
    }
}
```