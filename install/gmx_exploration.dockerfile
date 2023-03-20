# Run "docker build -f gmx_exploration.dockerfile --network=host --tag rid-gmx-exploration ." to build the image

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
    && conda install python=3.9 libtensorflow_cc=2.6.2=*cuda110* tensorflow=2.6.2=*cuda110* dpdata nccl -y

# Solve library inconsistency
RUN source activate base && rm ${CONDA_PREFIX}/lib/libtinfo.so* && ln -s /usr/lib/x86_64-linux-gnu/libtinfo.so.6 ${CONDA_PREFIX}/lib/libtinfo.so.6 \
    && rm ${CONDA_PREFIX}/lib/libcurl.so* && ln -s /usr/lib/x86_64-linux-gnu/libcurl.so.4 ${CONDA_PREFIX}/lib/libcurl.so.4

# Install necessary packages
RUN source activate base && pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple && pip install matplotlib==3.6.1 cython

# set enviroment variables
ENV LIBRARY_PATH "/opt/conda/lib:$LIBRARY_PATH"
ENV LD_LIBRARY_PATH "/opt/conda/lib:$LD_LIBRARY_PATH"

# Compile Plumed from source
RUN wget -O /root/plumed-2.8.1.tgz https://github.com/plumed/plumed2/releases/download/v2.8.1/plumed-2.8.1.tgz
RUN tar zxvf /root/plumed-2.8.1.tgz
#ADD plumed-2.8.1.tgz /root/
COPY DeePFE.cpp /root/plumed-2.8.1/src/bias/
RUN source activate base \
    && cd /root/plumed-2.8.1 \
    && ./configure --prefix=$CONDA_PREFIX \
                   CXXFLAGS="-O3 -I$CONDA_PREFIX/include/" \
                   LDFLAGS=" -L$CONDA_PREFIX/lib -ltensorflow_cc -ltensorflow_framework" \
    && make -j 6 \
    && make install \
    && rm -rf /root/plumed-2.8.1 /root/plumed-2.8.1.tgz

ENV PLUMED_KERNEL "/opt/conda/lib/libplumedKernel.@SOEXT@"

# Compile gromacs from source
RUN wget -O root/gromacs-2022.4.tar.gz https://github.com/gromacs/gromacs/archive/refs/tags/v2022.4.tar.gz --no-check-certificate
RUN tar zxvf /root/gromacs-2022.4.tar.gz
#ADD gromacs-2022.4.tar.gz /root/
RUN source activate base && echo $PLUMED_KERNEL && ls $CONDA_PREFIX/lib/ | grep plumed \
    && cd /root/gromacs-2022.4 \
    && echo -e "4\n" | plumed patch -p \
    && mkdir build \
    && cd build \
    && cmake .. -DGMX_BUILD_OWN_FFTW=ON \
                -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
                -DGMX_GPU=CUDA\
    && make -j 8 \
    && make install

#delete source files
RUN  rm -rf /root/gromacs-2022.4 /root/gromacs-2022.4.tar.gz \
    && rm -rf /var/lib/apt/lists/*

# Environment variables set manually since github actions will overwrite entrypoints
ENV CONDA_PREFIX=/opt/conda
RUN echo "source $CONDA_PREFIX/bin/GMXRC" >> /root/.bashrc

SHELL ["/bin/bash", "-c"]
