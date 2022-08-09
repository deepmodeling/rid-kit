# Configure minikube enviroment
**Start configuration**

Since rid-kit requires a large amount of resources, we recomment that you assign enough resources when staring the minikube enviroment. Otherwise there may be some unexpected errors after submitting the workflow. Depending on the machine resources, the recommended setting when starting the minikube is
```
minikube start --cpus 16 --disk-size=50g --extra-config=kubelet.max-pods=200 --image-mirror-country='cn'
```
the mirror image is for users in China. This assign 16 cpus, 50GB of disk and 100 max number of pods to the minikube node.

**You may encounter this error after starting the minikube**
```
    Unfortunately, an error has occurred:
            timed out waiting for the condition

    This error is likely caused by:
            - The kubelet is not running
            - The kubelet is unhealthy due to a misconfiguration of the node in some way (required cgroups disabled)

    If you are on a systemd-powered system, you can try to troubleshoot the error with the following commands:
            - 'systemctl status kubelet'
            - 'journalctl -xeu kubelet'

    Additionally, a control plane component may have crashed or exited when started by the container runtime.
    To troubleshoot, list all containers using your preferred container runtimes CLI.

    Here is one example how you may list all Kubernetes containers running in docker:
            - 'docker ps -a | grep kube | grep -v pause'
            Once you have found the failing container, you can inspect its logs with:
            - 'docker logs CONTAINERID'
```
This is expected, there might be an error related to pulling the image, specifically the "k8s.gcr.io/pause:3.6" image, the solution is as follows
```
minikube ssh
docker pull registry.aliyuncs.com/google_containers/pause:3.6
docker tag registry.aliyuncs.com/google_containers/pause:3.6 k8s.gcr.io/pause:3.6
exit
```
