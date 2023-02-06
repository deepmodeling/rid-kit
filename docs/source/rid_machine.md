# Setting Machine configuration of `RiD-kit`

`RiD-kit` uses a JSON-format file (typically `machine.json`) to configure machine resources to run the workflow. Here we explain different resources one by one.

`Rid-kit` currently supports three types of machine platform to do computation: `local machine`, `Bohrium cloud` and `Slurm machine`. `local machine` and `Bohrium cloud` use the ready-made `Docker image` to manage computation enviroment, while `Slurm machine` use `Conda` to manage computation enviroment.

## Local machine and Bohrium cloud
To run rid-kit workflow on a local machine or server, you need to build your own `k8s` enviroment by following the tutorials. But this is not recommended since `Rid-kit` workflow requires a lot of computation resources. 

To run rid-kit workflow on `Bohrium cloud`, you can either build your own `k8s` enviroment or simply use the community version of `k8s` [deepmodeling k8s](https://workflows.deepmodeling.com/).

Then the configuration is very easy, just configure different images to different rid-kit workflow steps. `Rid-kit` currently has six different `Docker images` to support different steps, add `registry.dp.tech/public/` if submitting to the `Bohrium cloud`.

`(registry.dp.tech/public/)pkufjhdocker/rid-gmx-exploration:stable` to do `run-exploration`. This image has gromacs compiled with plumed-patch, also supports tensorflow c++ interface, used to do neural network potential biased MD.

`(registry.dp.tech/public/)pkufjhdocker/rid-gmx-plumed:stable` to do `run-label`. This image has gromacs compiled with plumed-patch, used to do restrained or constrained MD.

`(registry.dp.tech/public/)pkufjhdocker/rid-gmx-tf:stable` to do `run-select` and `post-label`. This image has gromacs compiled with tensorflow python interface, used to select configuration in the `Exploration` step and calculating mean force in the constrained MD simulation. So if you are using restrained MD to calculate the mean force, `post-label` step can also use simple `pkufjhdocker/rid-tf-cpu:stable` images.

`(registry.dp.tech/public/)pkufjhdocker/rid-tf-gpu:stable` to do `train`. This image is the standard tensorflow image with GPU support, used to train the free-energy model.

`(registry.dp.tech/public/)pkufjhdocker/rid-tf-cpu:stable` to do other steps.

`(registry.dp.tech/public/)pkufjhdocker/gmx-exploration-dp:stable` to do `run-exploration` and `run-label`, if you want to run a DP potential molecular dynamics. This image has gromacs compiled with plumed-patch and DP potential support, also support tensorflow c++ interface, used to do neural network potential biased MD with DP potential.


### Local Machine Example
```JSON
{
    "resources": {
        "local_k8s_1": {
            "template_config" : {
                "image": "pkufjhdocker/rid-gmx-exploration:latest", 
                "image_pull_policy" : "IfNotPresent",
                "requests" : {"ephemeral-storage": "1Gi"}
            }
        },
        "local_k8s_2": {
            "template_config" : {
                "image": "pkufjhdocker/rid-gmx-plumed:latest", 
                "image_pull_policy" : "IfNotPresent",
                "requests" : {"ephemeral-storage": "1Gi"}
            }
        },
        "local_k8s_3": {
            "template_config" : {
                "image": "pkufjhdocker/rid-tf-cpu:latest", 
                "image_pull_policy" : "IfNotPresent",
                "requests" : {"ephemeral-storage": "1Gi"}
            }
        },
        "local_k8s_4": {
            "template_config" : {
                "image": "pkufjhdocker/rid-gmx-tf:latest", 
                "image_pull_policy" : "IfNotPresent",
                "requests" : {"ephemeral-storage": "1Gi"}
            }
        },
        "local_k8s_5": {
            "template_config" : {
                "image": "pkufjhdocker/rid-tf-gpu:latest", 
                "image_pull_policy" : "IfNotPresent",
                "requests" : {"ephemeral-storage": "1Gi"}
            }
        }
    },

    "tasks": {
        "prep_exploration_config": "local_k8s_3",
        "run_exploration_config": "local_k8s_1",
        "prep_label_config": "local_k8s_3",
        "run_label_config": "local_k8s_2",
        "post_label_config": "local_k8s_4",
        "prep_select_config": "local_k8s_3",
        "run_select_config": "local_k8s_4",
        "prep_data_config": "local_k8s_3",
        "run_train_config": "local_k8s_5",
        "workflow_steps_config": "local_k8s_3"
    }
}
```

### Bohrium example
```JSON
{
    "resources": {
        "bohrium1": {
            "template_config":{
                "requests" : {"ephemeral-storage": "1Gi"}
            },
            "executor":{
            "image": "dptechnology/dpdispatcher:latest",
            "machine_dict":{
                "batch_type": "DpCloudServer",
                "context_type": "DpCloudServerContext",
                "local_root" : "./",
                "remote_profile":{
                    "email": "",
                    "password": "",
                    "program_id": "",
                    "input_data":{
                        "api_version":2,
                        "job_type": "container",
                        "log_file": "tmp_log",
                        "job_name": "rid",
                        "scass_type":"c12_m92_1 * NVIDIA V100",
                        "platform": "ali",
                        "image_name":"registry.dp.tech/public/pkufjhdocker/rid-gmx-exploration:stable"
                        }
                }
            },
            "resources_dict":{
                "batch_type:": "Bohrium"
            }
        }
        },
        "bohrium2": {
            "template_config":{
                "requests" : {"ephemeral-storage": "1Gi"}
            },
            "executor":{
            "image": "dptechnology/dpdispatcher:latest",
            "machine_dict":{
                "batch_type": "DpCloudServer",
                "context_type": "DpCloudServerContext",
                "local_root" : "./",
                "remote_profile":{
                    "email": "",
                    "password": "",
                    "program_id": "",
                    "input_data":{
                        "api_version":2,
                        "job_type": "container",
                        "log_file": "tmp_log",
                        "job_name": "test_rid",
                        "scass_type":"c12_m92_1 * NVIDIA V100",
                        "platform": "ali",
                        "image_name":"registry.dp.tech/public/pkufjhdocker/rid-gmx-plumed:stable"
                        }
                }
            },
            "resources_dict":{
                "batch_type:": "Bohrium"
            }
        }
        },
        "bohrium3": {
            "template_config":{
                "requests" : {"ephemeral-storage": "1Gi"}
            },
            "executor":{
            "image": "dptechnology/dpdispatcher:latest",
            "machine_dict":{
                "batch_type": "DpCloudServer",
                "context_type": "DpCloudServerContext",
                "local_root" : "./",
                "remote_profile":{
                    "email": "",
                    "password": "",
                    "program_id": "",
                    "input_data":{
                        "api_version":2,
                        "job_type": "container",
                        "log_file": "tmp_log",
                        "job_name": "rid",
                        "scass_type":"c12_m92_1 * NVIDIA V100",
                        "platform": "ali",
                        "image_name":"registry.dp.tech/public/pkufjhdocker/rid-tf-cpu:stable"
                        }
                }
            },
            "resources_dict":{
                "batch_type:": "Bohrium"
            }
        }
        },
        "bohrium4": {
            "template_config":{
                "requests" : {"ephemeral-storage": "1Gi"}
            },
            "executor":{
            "image": "dptechnology/dpdispatcher:latest",
            "machine_dict":{
                "batch_type": "DpCloudServer",
                "context_type": "DpCloudServerContext",
                "local_root" : "./",
                "remote_profile":{
                    "email": "",
                    "password": "",
                    "program_id": "",
                    "input_data":{
                        "api_version":2,
                        "job_type": "container",
                        "log_file": "tmp_log",
                        "job_name": "rid",
                        "scass_type":"c12_m92_1 * NVIDIA V100",
                        "platform": "ali",
                        "image_name":"registry.dp.tech/public/pkufjhdocker/rid-gmx-tf:stable"
                        }
                }
            },
            "resources_dict":{
                "batch_type:": "Bohrium"
            }
        }
        },
        "bohrium5": {
            "template_config":{
                "requests" : {"ephemeral-storage": "1Gi"}
            },
            "executor":{
            "image": "dptechnology/dpdispatcher:latest",
            "machine_dict":{
                "batch_type": "DpCloudServer",
                "context_type": "DpCloudServerContext",
                "local_root" : "./",
                "remote_profile":{
                    "email": "",
                    "password": "",
                    "program_id": "",
                    "input_data":{
                        "api_version":2,
                        "job_type": "container",
                        "log_file": "tmp_log",
                        "job_name": "rid",
                        "scass_type":"c12_m92_1 * NVIDIA V100",
                        "platform": "ali",
                        "image_name":"registry.dp.tech/public/pkufjhdocker/rid-tf-gpu:stable"
                        }
                }
            },
            "resources_dict":{
                "batch_type:": "Bohrium"
            }
        }
        }
    },

    "tasks": {
        "prep_exploration_config": "bohrium3",
        "run_exploration_config": "bohrium1",
        "prep_label_config": "bohrium3",
        "run_label_config": "bohrium2",
        "post_label_config": "bohrium3",
        "prep_select_config": "bohrium3",
        "run_select_config": "bohrium4",
        "prep_data_config": "bohrium3",
        "run_train_config": "bohrium5",
        "workflow_steps_config": "bohrium3"
    }
}
```

## Slurm machine
To run rid-kit workflow on `Slurm machine`, you can either build your own `k8s` enviroment or simply use the community version of `k8s` [deepmodeling k8s](https://workflows.deepmodeling.com/).

If you want to submit `Rid-kit` workflow to `Slurm machine`, you have to compile the rid-kit enviroment yourself on the `Slurm machine`, we recommend to manage the enviroment in by `Conda`, you can check [Installation](install.md) for the compilation detail.

After compiling your `Conda` enviroment on the `Slurm machine`, the submission process is very easy now. Suppose the `Conda` enviroment can be activated by `rid.env` script, one submission example is the following:
### Slurm example
```JSON
{
    "resources": {
        "local_k8s": {
            "template_config" : {
                "image": "pkufjhdocker/rid-tf-cpu:latest", 
                "image_pull_policy" : "IfNotPresent",
                "requests" : {"ephemeral-storage": "1Gi"}
            }
        },
        "remote_slurm": {
            "template_config":{
                "requests" : {"ephemeral-storage": "1Gi"}
            },
            "executor":{
            "image": "dptechnology/dpdispatcher:latest",
            "machine_dict":{
                "batch_type": "Slurm",
                "context_type": "SSHContext",
                "local_root" : "./",
                "remote_root": "",
                "remote_profile":{
                    "hostname": "",
                    "username": "",
                    "password": "",
                    "port": "",
                    "timeout": 20
                }
            },
            "resources_dict":{
                "number_node": 1,
                "cpu_per_node": 8,
                "gpu_per_node": 1,
                "queue_name": "",
                "group_size": 1,
                "custom_flags": [
                    "#SBATCH --time=120:00:00",
                    "#SBATCH --exclude=3090-[000,007]"
                   ],
                "source_list": ["/path/to/rid.env"]
            }
        }
        }
    },

    "tasks": {
        "prep_exploration_config": "local_k8s",
        "run_exploration_config": "remote_slurm",
        "prep_label_config": "local_k8s",
        "run_label_config": "remote_slurm",
        "post_label_config": "remote_slurm",
        "prep_select_config": "local_k8s",
        "run_select_config": "remote_slurm",
        "prep_data_config": "local_k8s",
        "run_train_config": "remote_slurm",
        "workflow_steps_config": "local_k8s"
    }
}
```