# Installation of Software Environment

> ***Note***:
> 
> If you want to set environments by hand, please follow settings 1.1.1 ~ 1.1.4. 
You can also install environment of rid through the dockerfile, just run "docker build --net=host --tag rid-kit ." and then run "docker run -d -it --gpus all --net=host -v /mnt/workdir/:/mnt/workdir/ --shm-size=400g --name rid-container rid-kit" you can get the required enviroment.
## Install python and tensorflow

## Install tensorflow's C++ interface
### Using Bazel
It is relatively simple and one can follow the steps <https://docs.deepmodeling.com/projects/deepmd/en/master/install/install-tf.2.3.html>

The problem with Bazel compile is that it requires accessing website aboard. One way out is to use the compiled file and set path variables `PATH` and `LD_LIBRARY_PATH` accordingly.

### Using Conda
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


## Install plumed2.8.0
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

## Install gromacs 2021.4

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
