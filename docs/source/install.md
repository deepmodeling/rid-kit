# Environment settings

The `enviroment` of rid-kit software is a bit complex, it uses `dflow` to manage the workflow in `kubernetes` enviroment refered to as `management environment`. The actual rid-kit code runs on local or remote (slurm, cloud) enviroment refered to as `computation environment`. The settings of the `management enviroment` has already been described in [Tutorial](tutorial.ipynb), only settings of `computation enviroment` is given as follows.

## Configuration of computation environment
> ***Note***:
> 
> If you want to set environments by hand, please follow settings of this section. 
>
> You can also set computation environment of rid through docker, just run "docker pull pkufjhdocker/rid-kit:latest" to get the docker image we built. Or build the image through the dockerfile in rid-kit/install directory.

### Install tensorflow's C++ interface and other necessary packages
We recommend using conda to manage the python enviroment. 
Use the command
```
conda create -n rid-drop python=3.9 libtensorflow_cc=2.6.2=*cuda110* tensorflow=2.6.2=*cuda110* nccl mdtraj numpy
```
to get the compiled libtensorflow_cc and other necessary packages.
After installation, activate the enviroment
```
conda activate rid-drop
```

### Install plumed2.8.0
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
export LIBRARY_PATH=${CONDA_PREFIX}/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:$LD_LIBRARY_PATH
export PLUMED_KERNEL=${CONDA_PREFIX}/lib/libplumedKernel.@SOEXT@
```

### Install gromacs 2021.4

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
source ${CONDA_PREFIX}/bin/GMXRC
```
