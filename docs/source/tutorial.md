# Tutorial

In this tutorial, we will learn how to deploy `RiD-kit` with your own `Kubenetes` enviroment and run a simple case of alanine dipeptide.


## Installation of `dflow` and `rid-kit`
With the power of `dflow`, users can easily minitor the whole workflow of RiD tasks and dispatch their tasks to various computational resources. Before you use it, you should have `dflow` installed on your host computer (your PC or a remote server). 

It it necessary to emphasize that, the computational nodes and monitor nodes are seperated. With `dflow`, you can deploy `dflow` and `rid` on your PC and achieve expensive computation on other resources (like `Slurm` and Cloud Platform) without any further effort.

Instructions of `dflow` installation are provided in detail on its [Github page](https://github.com/deepmodeling/dflow#Installdflow). Prerequisites of `dflow` usage are `Docker` and `Kubenetes`, where their main pages ([Docker](https://docs.docker.com/engine/install/) & [Kubenetes](https://kubernetes.io/docs/tasks/tools/)) include how you can install them. Besides, `dflow` repo also provides with easy-install shell scripts on [dflow/scripts](https://github.com/deepmodeling/dflow/tree/master/scripts) to install `Docker` & `Kubenetes` & `dflow` and make port-forwarding.

Instead of using the easy-install scripts, it is recommended to configure your own minikube enviroment

Installation of `minikube` is very easy, just follow the page [Minikube](https://minikube.sigs.k8s.io/docs/start/) based on your machine archetecture.

Then use minikube to start a k8s environment. This can be done by just type `minikube start`, but since running rid-kit requires relatively huge amount of resources, it is recommended to allocate more resources when executing `minikube start`, just as follows.

> **Note:**
>
> Don't try to run `minikube` with root privileges, otherwise an error may occur:
>
> `Exiting due to DRV_AS_ROOT: The "docker" driver should not be used with root privileges.`
>

```bash
#choose the location of the minikube enviroment
export MINIKUBE_HOME=~/.minikube  
#allocate enough memory based on your machine in case of high parallelism
minikube start --cpus 8 --memory 8192mb --kubernetes-version=1.23.9 --image-mirror-country='cn'
# mount the storage path if you are using machine with shared memory storage system, change the minio host path accordingly
minikube start --cpus 8 --memory 8192mb --kubernetes-version=1.23.9 --mount --mount-string="/path_on_your_machine:/data2" --image-mirror-country='cn'
```

A further step to  configure `argo` service is to run:

```bash
kubectl create ns argo
kubectl apply -n argo -f https://raw.githubusercontent.com/deepmodeling/dflow/master/manifests/quick-start-postgres.yaml
```

If you want to mount your local path into the minio path, change the minio host path in the yaml file

```yaml
      - hostPath:
          path: /data2/minio
          type: DirectoryOrCreate
```

Now you should have installed `Docker` and `minikube` properly. Run commands to check their status. For `minikube`, you should wait util all servers keep `running`. This may take a couple of minutes.

```bash
# check minikube status
minikube status
# check k8s status in argo namespace
kubectl get pods -n argo
```

### Installation of RiD-kit

Now we install `rid-kit` on the host machine. To meet the minimum requirments, the prerequisites of third-party python package should be installed:

* tensorflow-cpu or gpu
* mdtraj
* numpy
* scikit-learn
* pydflow
* google
* mdtraj
* dpdata
* parmed
* dpdispatcher
* lbg
* matplotlib

which are also listed in `rid-kit/requirements.txt`. Then change directory to `rid-kit` repo and run:

```bash
# the rid-kit repo path
git checkout dflow
pip install setuptools_scm
pip install .
```

## Configuration of Computational Environment

In RiD workflow, `dflow` helps send computation tasks to resources with peoper environment configured. 

There are four main modules and several workflow steps in RiD procedures and each module or step needs different environments:

* `Exploration/Sampling`:  `Gromacs`, `PLUMED2` modified by `DeepFE.cpp`, `Tensorflow` C++ interface. (prefer GPU)
* `Selection`:  `Tensroflow` Python interface.
* `Labeling`:  `Gromacs`, `PLUMED2`. (prefer GPU)
* `Training`:  `Tensroflow` Python interface. (prefer GPU).
* `Workflow steps`:  `Python`.

`dflow` supports different resources including `Slurm` clusters, `K8S` local machines and `Cloud Server`.

* For `Slurm`, configure computational environments on your `Slurm` following the instructions in [Environment settings](https://github.com/PKUfjh/rid-kit/blob/dflow/docs/source/install.md). With `dflow`, `rid-kit` send tasks to `Slurm` nodes from the host machines remotely without manually logging in the cluster.
* For local resources, just use the docker images we have built. No further manual configuration needed. We also provide `Dockerfile` of our images to enable flexible modification.
* For `Cloud Server`, like `Bohrium`, use public images and no further manual configuration needed.

We highly recommend using `Bohrium` to do computation since no enviroment installation is required, you only need to use the public image.Specifically, the public images for rid-kit involves the following: `registry.dp.tech/public/pkufjhdocker/rid-gmx-exploration:stable` for doing `Exploration`, `registry.dp.tech/public/pkufjhdocker/rid-gmx-plumed:stable` for doing `Labeling`, `registry.dp.tech/public/pkufjhdocker/rid-gmx-tf:stable` for doing `Selection`, `registry.dp.tech/public/pkufjhdocker/rid-tf-gpu:stable` for `Training` and `registry.dp.tech/public/pkufjhdocker/rid-tf-cpu:stable` for doing other steps.

## Prepare machine configuration JSON.

`rid-kit` uses `JSON` file to manage resources. In `machine.json`, define your own `resources` and dispatch `tasks` to them.
* In key `resources`, you define your own resources types. Resource names and their numbers are custom.
* In key `tasks`, you distribute resources you have defined to tasks of RiD. Do not change task names in `tasks` as they are fixed in codes.

Generally, we would like to run low-cost tasks on cpu nodes or locally and submit high-cost tasks to `Slurm` or `Clouds`. 

If submitting to mixed `Slurm` and `local` environment, a `machine.json` example is in [rid/template/machine_slurm_k8s.json](https://github.com/PKUfjh/rid-kit/blob/dflow/rid/template/machine_slurm_k8s.json):

If submitting to `Bohrium` cloud platform, a `machine.json` example is in [rid/template/machine_bohrium.json](https://github.com/PKUfjh/rid-kit/blob/dflow/rid/template/machine_bohrium.json).

If submitting to mixed `Bohrium` and `local` platform, a `machine.json` example is in [rid/template/machine_bohrium_k8s.json](https://github.com/PKUfjh/rid-kit/blob/dflow/rid/template/machine_bohrium_k8s.json).

If submitting to `Slurm` environment without `k8s` enviroment, a `machine.json` example is in [rid/template/machine_slurm_local.json](https://github.com/PKUfjh/rid-kit/blob/dflow/rid/template/machine_slurm_local.json):

If you submit jobs to `Slurm` enviroment, you will have to compile the computation enviroment on `Slurm` server for yourself. While in `Bohrium` cloud platform the computation images are all set, no extra compilation is needed. So we recommend to submit jobs to `Cloud` enviroment like `Bohrium`.

For more details about setting the machine resources, check [Rid machine resourse](https://github.com/PKUfjh/rid-kit/blob/dflow/docs/source/rid_machine.md)

## Prepare rid configuration JSON.
If you want to use `gromacs` as sampler, [rid/template/rid_gmx_dih.json](https://github.com/PKUfjh/rid-kit/blob/dflow/rid/template/rid_gmx_dih.json) is an example.

If you want to use `gromacs` with `DeepMD` potential as sampler, [rid/template/rid_gmx_dp.json](https://github.com/PKUfjh/rid-kit/blob/dflow/rid/template/rid_gmx_dp.json) is an example.

If you want to use `lammps` with `DeepMD` potential as sampler, [rid/template/rid_lmp_dp.json](https://github.com/PKUfjh/rid-kit/blob/dflow/rid/template/rid_lmp_dp.json) is an example.

For more details about setting the rid configuration, check [Rid configuration](https://github.com/PKUfjh/rid-kit/blob/dflow/docs/source/rid_configuration.md)

## Get Started!

Assume you have learn the basic knowledge of reinforced dynamics which we won't describe again here. 

Users can monitor workflows from browser UI. To enable that, you should forward ports of `argo` and `minio`. These could be achieved by `rid port-forward`.

```bash
rid port-forward
```

In this case, we try to explore the phase space of alanine dipeptide. Prepare your initial conformation files in `.gro` format, topology file in `.top` format and configuration file `rid.json`. For convenience, we have prepared on at `rid/template/rid.json`. 
Remember also provide your own forcefield files. Collect all these files into a directory and feed its path to `rid-kit` by flag `-i`.

A minimum case was prepared in `rid-kit/tests/data/000`. After configurating your `machine.json`, then run `rid submit`:

```bash
rid submit -i ./tests/data/000 -c ./rid/template/rid_gmx_dih.json -m ./rid/template/machine_bohrium_k8s.json
```

You can specify the workflow name by providing WORKFLOW_ID after "-d", for example:

```bash
rid submit -i ./tests/data/000 -c ./rid/template/rid_gmx_dih.json -m ./rid/template/machine_bohrium_k8s.json -d ala-dipeptide-1
```

**Note that the defined workflow-id should only contain lower case alphanumeric character, and specifal character "-".**

INFO indicates that this task has been submitted succussfully. Record this workflow ID as we may use it later.

Visit the `url` given by the last two lines, all workflows and corresponding files are listed on UI.

Command lines are also supported. Run `rid ls` to list your workflows and their status.

```bash
! rid ls
```

`rid-kit` is based on `dflow`, `argo` and `minikube`. So further complex and flexible managements of workflows can be achieved by their command lines. like `kubectl get pods -n argo` and `argo show`.

For failed tasks, you may want to remove them or resubmit them from the failure steps.

For `remove`:

```bash
# rid rm task-ID
rid rm reinforced-dynamics-bsc7j 
```

For `resubmit` to modify and continue workflow:

```bash
# suppose the original workflow id is OLD_ID
rid resubmit -i your_dir -c path_to_rid.json -m path_to_machine.json OLD_ID -d NEW_ID
```

If you want to resubmit from a particular `iteration` and `step`:

```bash
rid resubmit -i your_dir -c path_to_rid.json -m path_to_machine.json OLD_ID -t ITERATION-ID -p STEP-KEY -d NEW_ID
```

The `Workflow-ID` is something like `reinforced-dynamics-jq4jn` appeared in the argo UI. If the workflow is archived, its name will appear as somethinig like `a8463748-e15a-4f2c-882b-bfd981a76dac`. If you want to resubmit a archived workflow, you have to provide the archived ID rather than the initial ID. The `ITERATION-ID` is just `n`th iteration the workflow has been executed. The `STEP-KEY` in rid includes the following steps: `prep-exploration`, `run-exploration`, `prep-select`, `run-select`, `prep-label`, `run-label`, `label-stats`, `collect-data`, `merge-data`, `train`, `model-devi`.
