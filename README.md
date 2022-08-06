
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
	* 2.3. [Install rid package](#Installridpackage)
* 3. [Quick Start](#QuickStart)
* 4. [Main procedure of RiD](#MainprocedureofRiD)
		* 4.1. [a. Biased MD](#a.BiasedMD)
		* 4.2. [b. Restrained MD](#b.RestrainedMD)
		* 4.3. [c. Neural network training](#c.Neuralnetworktraining)


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

### 2.3. <a name='Installridpackage'></a>**Install rid package**
Now you have all dependence of RiD (Gromacs, Tensorflow, and a conda environment).
~~~bash
cd rit-kit
python setup.py install
~~~


##  3. <a name='QuickStart'></a>**Quick Start**

Please see `example/tutorial.ipynb` for detail.


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