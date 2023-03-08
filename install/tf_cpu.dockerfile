
# Build from tensorflow/tensorflow
FROM tensorflow/tensorflow:latest
SHELL ["/bin/bash", "-c"]

RUN pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple \
    && pip install dpdata matplotlib mdtraj scikit-learn parmed

RUN apt-get install -y python3-tk
SHELL ["/bin/bash", "-c"]