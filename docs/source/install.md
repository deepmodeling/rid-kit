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
conda create -n rid-drop python=3.9 libtensorflow_cc=2.6.2=*cuda110* tensorflow=2.6.2=*cuda110* nccl mdtraj matplotlib==3.6.1 scikit-learn cmake dpdata parmed cython -c conda-forge
```
to get the compiled libtensorflow_cc and other necessary packages.
After installation, activate the enviroment and set library path
```bash
conda activate rid-drop
export LIBRARY_PATH=${CONDA_PREFIX}/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:$LD_LIBRARY_PATH
```

### Install plumed2.8.1
To compile plumed with our free energy model interface, you need to copy compiled `DeePFE.cpp` to the plumed directory. This file locates at `rid-kit/install/DeePFE.cpp`.
Note that here we use c++ 14 standard to compile, the gcc version need to be newer than 6.1.0. In slurm enviroment, if the software is managed by `Enviroment Module`, one can load gcc such as `module load gcc/10.1.0` or other available versions.

set CC and CXX enviroment variables
```
export CC=/yourpath/to/gcc-10.1.0
export CXX=/yourpath/to/g++-10.1.0
```

```bash
tar -xvzf plumed-2.8.1.tgz
cp DeePFE.cpp plumed-2.8.1/src/bias
cd plumed-2.8.1
./configure --prefix=$CONDA_PREFIX \
                   CXXFLAGS="-O3 -I$CONDA_PREFIX/include/" \
                   LDFLAGS=" -L$CONDA_PREFIX/lib -ltensorflow_cc -ltensorflow_framework"
make -j 6
make install
```
Set the enviroment variables
```bash
export PLUMED_KERNEL=${CONDA_PREFIX}/lib/libplumedKernel.@SOEXT@
```

### Install gromacs 2022.4
Note that in order to compile gromacs with GPU support, you need to have CUDA toolkit installed. In slurm enviroment, if the software is managed by `Enviroment Module`, one can load CUDA toolkits such as `module load cuda/11.4`. Also note that in most slurm environment, there is no network connection at computing node, but if you want to compile gromacs without fftw installed, you will have to compile it on login-in node where network connection is on.

However if the network is not on in the slurm machine, you can first compile fftw manually by
```bash
cd fftw-3.3.8
./configure --prefix=$CONDA_PREFIX
make
make install
```

```bash
tar -xzvf gromacs-2022.4.tar.gz
cd gromacs-2022.4
plumed patch -p # make sure you use the correct gromacs version
mkdir build
cd build
# if the network is on
cmake .. -DGMX_BUILD_OWN_FFTW=ON \
                -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
                -DGMX_GPU=CUDA
# if the network is off
cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
            -DGMX_GPU=CUDA
make -j 4
make install
```
Set the enviroment variables
```bash
source ${CONDA_PREFIX}/bin/GMXRC
```

### Install Lammps-stable_23Jun2022_update1
```bash
tar -xzvf stable_23Jun2022_update1.tar.gz
cd lammps-stable_23Jun2022_update1/src
```
make packages you will use
```bash
make yes-extra-fix
make yes-extra-dump
make yes-kspace
```
make plumed
```bash
make lib-plumed args="-p $CONDA_PREFIX -m runtime"
make yes-plumed
make serial -j4 #build serial version of lammps
ln -s /home/dongdong/software/lammps-stable_23Jun2022_update1/src/lmp_serial /home/dongdong/software/anaconda3/envs/rid_lmp/bin/lmp_serial
```
