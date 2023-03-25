
# Build from tensorflow/tensorflow gpu version
FROM tensorflow/tensorflow:devel-gpu
SHELL ["/bin/bash", "-c"]

RUN pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple \
    && pip install dpdata matplotlib mdtraj scikit-learn

RUN apt-get install -y python3-tk
SHELL ["/bin/bash", "-c"]