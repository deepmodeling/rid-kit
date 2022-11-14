FROM dptechnology/dflow:latest

WORKDIR /data/rid-kit
COPY ./ ./
RUN pip config set global.index-url http://mirrors.aliyun.com/pypi/simple
RUN pip config set install.trusted-host mirrors.aliyun.com
RUN pip install setuptools_scm
RUN pip install pyparsing
RUN pip install .
RUN pip install --upgrade protobuf==3.20.1