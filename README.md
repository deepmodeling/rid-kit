# Table of contents
- [About Rid-kit](#about-rid-kit)
- [Main procedure of RiD](#main-procedure-of-rid)
- [Use Rid-kit](#use-rid-kit)
- [Installation of enviroment](#installation-of-enviroment)
- [Installation of DeepMD potential](#installation-of-dp)
- [Configure simulations](#configure-simulations)
- [Troubleshooting](#troubleshooting)

# About Rid-kit
Rid-kit is a package written in Python, designed to do enhanced sampling using reinforced dynamics. It aims to learn the free energy surface on the fly during MD run, and uses it as the bias potential during the next MD run. Its advantage is the ability to use a large number of CVs (100), thus can be used to simulate conformational changes of big molecules such as in the problem of protein folding.

For more information, check the [documentation](https://rid-kit-dflow.readthedocs.io/en/latest/).

# Main procedure of RiD

RiD will run in iterations. Every iteration contains tasks below:

1. Biased MD;
2. Restrained MD;
3. Training neural network.

## Biased MD

Just like Metadynamics, RiD will sample based on a bias potential given by NN models. An uncertainty indicator will direct the process of adding bias potential.

## Restrained MD

This procedure will calculate the mean force based on the sampling results, which can generate data set for training. 

## Neural network training

A fully connected NN will be trained via sampling data. This network will generate a map from selected CV to free energy.

A more detailed description of RiD is published now, please see:

>  [1]  Zhang, L., Wang, H., E, W.. Reinforced dynamics for enhanced sampling in large atomic and molecular systems[J]. The Journal of chemical physics, 2018, 148(12): 124113.
>  
>  [2]  Wang, D., Wang, Y., Chang, J., Zhang, L., & Wang, H. (2022). Efficient sampling of high-dimensional free energy landscapes using adaptive reinforced dynamics. Nature Computational Science, 2(1), 20-29.

# Use Rid-kit

A quick-start on using Rid-kit can be found as follows:

- [Tutorial](/docs/source/tutorial.ipynb)

# Enviroment settings

Installation of the computation enviroment can be found in
- [Installation](docs/source/install.md)

# Installation-of-dp
Installation of the DeepMD potential support can be found in
- [Installation of dp](docs/source/install_dp.md)
  
# Configure simulations
- [configure simulations](docs/source/rid_configuration.md)

# Workflow Synopsis

- ![image](docs/pics/rid_workflow.jpg)

# Troubleshooting
- [Installation](docs/source/troubleshooting/installation.md)
- [Minikube](docs/source/troubleshooting/minikube.md)
- [Submit](docs/source/troubleshooting/submit.md)
