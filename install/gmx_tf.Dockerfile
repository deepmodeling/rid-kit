# Run "docker build --network=host --tag *** ." to build the image

# Build from gromacs/gromacs iamge
FROM gromacs/gromacs:2022.2
SHELL ["/bin/bash", "-c"]

# Add public key of the image to trusted set of keys
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys A4B469963BF863CC

# Change mirrors for apt-get to aliyun 
RUN sed -i s@/archive.ubuntu.com/@/mirrors.aliyun.com/@g /etc/apt/sources.list

# Install essentials
RUN DEBIAN_FRONTEND=noninteractive apt-get update\
    && DEBIAN_FRONTEND=noninteractive apt-get install -y pip python3-tk

# Install python packages
RUN pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple \
    && pip install tensorflow dpdata mdtraj scikit-learn numpy parmed

SHELL ["/bin/bash", "-c"]