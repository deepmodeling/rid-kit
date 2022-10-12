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
conda create -n rid-dp python=3.9 libtensorflow_cc=2.6.2=*cuda110* tensorflow=2.6.2=*cuda110* nccl mdtraj numpy scikit-learn cmake -c conda-forge
```
to get the compiled libtensorflow_cc and other necessary packages.
After installation, activate the enviroment and set library path
```bash
conda activate rid-dp
export LIBRARY_PATH=${CONDA_PREFIX}/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:$LD_LIBRARY_PATH
```

## Install the DeepMD-kit's C++ interface
Now go to the source code directory of DeePMD-kit and make a build place.
```
cd $deepmd_source_dir/source
mkdir build 
cd build
```

install DeepMD into conda path
```
cmake -DTENSORFLOW_ROOT=$CONDA_PREFIX -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DUSE_CUDA_TOOLKIT=TRUE ..
```
You should use CUDA version>= 11.0, otherwise the following error will occur
```
/home/software/deepmd-kit/source/lib/src/cuda/neighbor_list.cu:4:10: 致命错误：cub/block/block_scan.cuh：没有那个文件或目录
 #include <cub/block/block_scan.cuh>
          ^~~~~~~~~~~~~~~~~~~~~~~~~~
编译中断。
CMake Error at deepmd_op_cuda_generated_neighbor_list.cu.o.release.cmake:222 (message):
  Error generating
  /home/software/deepmd-kit/source/build/lib/src/cuda/CMakeFiles/deepmd_op_cuda.dir//./deepmd_op_cuda_generated_neighbor_list.cu.o
```

In slurm enviroment, if the software is managed by `Enviroment Module`, one can load CUDA toolkits such as `module load cuda/11.1`.
If the cmake has been executed successfully, then run the following make commands to build the package:
```
make -j4
make install
```
If everything works fine, you should find executable `dp_gmx_patch` in $CONDA_PREFIX.

### Install plumed2.6.1
You need to copy compiled `DeePFE.cpp` to the plumed directory. This file locates at `rid-kit/install/DeePFE.cpp`.
Note that here we use c++ 14 standard to compile, the gcc version need to be newer than 6.1.0. In slurm enviroment, if the software is managed by `Enviroment Module`, one can load gcc such as `module load gcc/7.2.0` or other available versions.

set CC and CXX enviroment variables
```
export CC=/yourpath/to/gcc-7.2.0
export CXX=/yourpath/to/g++-7.2.0
```

```bash
tar -xvzf plumed-2.6.1.tgz
cp DeePFE.cpp plumed-2.6.1/src/bias
cd plumed-2.6.1
./configure --prefix=$CONDA_PREFIX \
                   --enable-cxx=14 \
                   CXXFLAGS="-std=gnu++14 -I $CONDA_PREFIX/include/" \
                   LDFLAGS=" -L$CONDA_PREFIX/lib -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=$CONDA_PREFIX/lib/" \
make -j 6
make install
```
Set the enviroment variables
```bash
export PLUMED_KERNEL=${CONDA_PREFIX}/lib/libplumedKernel.@SOEXT@
```

### Install gromacs 2020.2
Note that in order to compile gromacs with GPU support, you need to have CUDA toolkit installed. In slurm enviroment, if the software is managed by `Enviroment Module`, one can load CUDA toolkits such as `module load cuda/11.1`. Also note that in most slurm environment, there is no network connection at computing node, but if you want to compile gromacs without fftw installed, you will have to compile it on login-in node where network connection is on.

```bash
tar -xzvf gromacs-2020.2.tar.gz
cd gromacs-2020.2
plumed patch -p  # make sure you use the correct plumed version
dp_gmx_patch -d $gromacs_root -v $version -p #version=2020.2 is supported
mkdir build
cd build
export CMAKE_PREFIX_PATH="/path/to/fftw-3.3.9" # fftw libraries
cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
                -DGMX_GPU=CUDA  \
                -DGMX_SIMD=avx_512
make -j 4
make install
```
Set the enviroment variables
```bash
source ${CONDA_PREFIX}/bin/GMXRC
```
