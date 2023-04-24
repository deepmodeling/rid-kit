# Run "docker build -f lmp_dp.dockerfile --network=host --tag *** ." to build the image

# Build from nvidia/cuda iamge
FROM nvidia/cuda:11.0.3-cudnn8-devel-ubuntu20.04
SHELL ["/bin/bash", "-c"]

# Add public key of the image to trusted set of keys
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys A4B469963BF863CC

# Change mirrors for apt-get to aliyun 
RUN sed -i s@/archive.ubuntu.com/@/mirrors.aliyun.com/@g /etc/apt/sources.list

# Install essentials
RUN DEBIAN_FRONTEND=noninteractive apt-get update\
    && DEBIAN_FRONTEND=noninteractive apt-get install -y software-properties-common \
    && DEBIAN_FRONTEND=noninteractive add-apt-repository -y ppa:ubuntu-toolchain-r/test \
    && DEBIAN_FRONTEND=noninteractive apt-get update \
    && echo -e "6\n1\n" | DEBIAN_FRONTEND=noninteractive apt-get install -y wget gcc g++ cmake vim libstdc++6 python3-tk
ENV CC=/usr/bin/gcc
ENV CXX=/usr/bin/g++

# Install Miniconda package manager.
#RUN wget -O /root/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
ADD Miniconda3-latest-Linux-x86_64.sh /root/miniconda.sh
RUN bash /root/miniconda.sh -b -p /opt/conda
RUN rm /root/miniconda.sh
ENV PATH /opt/conda/bin:$PATH
RUN conda init bash && source /root/.bashrc \
    && conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/  \
    && conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/  \
    && conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch/  \
    && conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge  \
    && conda install python=3.9 libtensorflow_cc=2.6.2=*cuda110* tensorflow=2.6.2=*cuda110* dpdata nccl numpy -y

# Solve library inconsistency
RUN source activate base && rm ${CONDA_PREFIX}/lib/libtinfo.so* && ln -s /usr/lib/x86_64-linux-gnu/libtinfo.so.6 ${CONDA_PREFIX}/lib/libtinfo.so.6 \
    && rm ${CONDA_PREFIX}/lib/libcurl.so* && ln -s /usr/lib/x86_64-linux-gnu/libcurl.so.4 ${CONDA_PREFIX}/lib/libcurl.so.4

ENV LIBRARY_PATH "/opt/conda/lib:$LIBRARY_PATH"
ENV LD_LIBRARY_PATH "/opt/conda/lib:$LD_LIBRARY_PATH"

# Compile deepmd-kit
ADD ./deepmd-kit /root/deepmd-kit
RUN cd /root/deepmd-kit/source\
    && mkdir ./build\
    && cd ./build\
    && cmake -DTENSORFLOW_ROOT=/opt/conda -DCMAKE_INSTALL_PREFIX=/opt/conda -DUSE_CUDA_TOOLKIT=TRUE ..\
    && make -j4\
    && make install\
    && make lammps

# Compile Plumed from source
# RUN wget -O /root/plumed-2.8.0.tgz https://github.com/plumed/plumed2/releases/download/v2.8.0/plumed-2.8.0.tgz
# RUN tar zxvf /root/plumed-2.8.0.tgz
ADD plumed-2.8.0.tgz /root/
COPY DeePFE.cpp /root/plumed-2.8.0/src/bias/
RUN source activate base \
    && cd /root/plumed-2.8.0 \
    && ./configure --prefix=$CONDA_PREFIX \
                   --enable-cxx=14 \
                   CXXFLAGS="-std=gnu++14 -I $CONDA_PREFIX/include/" \
                   LDFLAGS=" -L$CONDA_PREFIX/lib -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=$CONDA_PREFIX/lib/" \
    && make -j 6 \
    && make install \
    && rm -rf /root/plumed-2.8.0 /root/plumed-2.8.0.tgz

ENV PLUMED_KERNEL "/opt/conda/lib/libplumedKernel.so"

# Compile gromacs from source
# RUN wget -O root/gromacs-2021.4.tar.gz https://github.com/gromacs/gromacs/archive/refs/tags/v2021.4.tar.gz --no-check-certificate
# RUN tar zxvf /root/gromacs-2021.4.tar.gz
ADD stable_23Jun2022_update1.tar.gz /root/
RUN source activate base && echo $PLUMED_KERNEL && ls $CONDA_PREFIX/lib/ | grep plumed \
    && cd /root/lammps-stable_23Jun2022_update1/src\
    && cp -r /root/deepmd-kit/source/build/USER-DEEPMD .\
    && make yes-user-deepmd\
    && make yes-extra-fix\
    && make yes-extra-dump\
    && make yes-kspace\
    && make yes-molecule\
    && make yes-rigid\
    && make lib-plumed args="-p $CONDA_PREFIX -m runtime"\
    && make yes-plumed\
    && make serial -j4

# manage the executables and source code
RUN source activate base && ln -s /root/lammps-stable_23Jun2022_update1/src/lmp_serial $CONDA_PREFIX/bin/lmp_serial

SHELL ["/bin/bash", "-c"]