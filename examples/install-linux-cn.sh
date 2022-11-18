#!/bin/bash

function INFO() {
    echo "[INFO] $@"
}

function WARNING() {
    echo >&2 "[WARNING] $@"
}

function ERROR() {
    echo >&2 "[ERROR] $@"
}

docker_path=$(which docker)
if [[ -n "$docker_path" ]]; then
    INFO "Found docker executable at $docker_path"
else
    INFO "Docker not found, installing docker..."
    sh -c "$(curl -fsSL https://get.docker.com/)"
    if [[ $? != 0 ]]; then
        ERROR "Fail to install docker"
        exit 1
    fi
fi

minikube_path=$(which minikube)
if [[ -n "$minikube_path" ]]; then
    INFO "Found minikube binary at $minikube_path"
else
    INFO "Minikube not found, installing minikube 1.25.2 ..."
    curl -o minikube -L https://registry.npmmirror.com/-/binary/minikube/v1.25.2/minikube-linux-amd64
    if [[ $? != 0 ]]; then
        ERROR "Fail to download minikube"
        exit 1
    fi
    sudo install minikube /usr/local/bin/minikube
    if [[ $? != 0 ]]; then
        ERROR "Fail to install minikube"
        exit 1
    fi
fi

kubectl=$(which kubectl)
if [[ $? != 0 ]]; then
    curl -LO "https://dl.k8s.io/release/$(curl -L -s https://dl.k8s.io/release/stable.txt)/bin/linux/amd64/kubectl"
    if [[ $? != 0 ]]; then
        ERROR "Fail to download kubectl"
        exit 1
    fi
    sudo install kubectl /usr/local/bin/kubectl
    if [[ $? != 0 ]]; then
        ERROR "Fail to install kubectl"
        exit 1
    fi
fi


# if [[ $EUID > 0 ]]; then
#     minikube start --image-mirror-country=cn $@
# else 
#     INFO "minikube can not start with root user"
#     INFO "Creating new user"
#     adduser developer
#     usermod -aG sudo developer
#     INFO "Login to the newly created user"
#     su - developer
#     INFO "Add user to the Docker Group"
#     sudo groupadd docker
#     sudo usermod -aG docker $USER
#     exit
#     su - developer
# fi 




minikube status 1>/dev/null 2>/dev/null
if [[ $? < 2 ]]; then
    INFO "Minikube has been started"
else
    INFO "Starting minikube..."
    minikube start --image-mirror-country=cn $@
    if [[ $? != 0 ]]; then
        ERROR "Fail to start minikube"
        exit 1
    fi
fi

kubectl create ns argo 1>/dev/null 2>/dev/null
wget https://raw.githubusercontent.com/deepmodeling/dflow/master/manifests/quick-start-postgres-stable-cn.yaml
kubectl apply -n argo -f quick-start-postgres-stable-cn.yaml 1>/dev/null
if [[ $? != 0 ]]; then
    ERROR "Fail to apply argo yaml"
    exit 1
fi

function waitForReady() {
    while true; do
        ready=$(kubectl get deployment $1 -n argo -o jsonpath='{.status.readyReplicas}')
        replicas=$(kubectl get deployment $1 -n argo -o jsonpath='{.status.replicas}')
        if [[ $? != 0 ]]; then
            ERROR "Fail to get status of $1"
            exit 1
        fi
        if [[ $replicas > 0 && $ready == $replicas ]]; then
            INFO "$1 has been ready..."
            break
        fi
        INFO "Waiting for $1 ready..."
        sleep 3
    done
}

waitForReady argo-server
waitForReady httpbin
waitForReady minio
waitForReady postgres
waitForReady workflow-controller

function forward() {
    pid=`ps -ef | grep port-forward | grep $1 | grep $2 | awk '{print $2}'`
    if [[ -n "$pid" ]]; then
        kill -9 $pid
    fi
    INFO "Forwarding $1:$2 to localhost:$2"
    nohup kubectl -n argo port-forward deployment/$1 $2:$2 --address 0.0.0.0 &
}
forward argo-server 2746
forward minio 9000
forward minio 9001

sleep 3
INFO "dflow server has been installed successfully!"
