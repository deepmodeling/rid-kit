# Submit
**You may encounter the error when submit the workflow**
```
...
urllib3.exceptions.MaxRetryError: HTTPConnectionPool(host='127.0.0.1', port=9000): Max retries exceeded with url: /my-bucket?location= (Caused by NewConnectionError('<urllib3.connection.HTTPConnection object at 0x7fdab9a05760>: Failed to establish a new connection: [Errno 111] Connection refused'))
```
The reason for this error is that the port has not been forwarded in a proper way, just type
```
rid port-forward
```
If you see something like this
```
2022-08-09 14:12:07 | INFO | rid.entrypoint.server | Port "agro-server" has been launched and running.
2022-08-09 14:12:07 | INFO | rid.entrypoint.server | Port "minio-server" has been launched and running.
2022-08-09 14:12:07 | INFO | rid.entrypoint.server | Port "minio-ui" has been launched and running.
```
then you are done.

**You may encounter this message on argo UI after submitting the workflow to Slurm**
```
ImagePullBackOff: Back-off pulling image "dptechnology/dflow-extender"
```
The status "ImagePullBackOff" means that a Pod couldn’t start, because Kubernetes couldn’t pull a container image. The ‘BackOff’ part means that Kubernetes will keep trying to pull the image, with an increasing delay (‘back-off’). When this problem happens, you can pull the image on your host machine (which started your minikube environment), then you can type
```
minikube image load dptechnology/dflow-extender
```
to load the image to the minikube enviroment. Furthermore, since the default setting of image-pull policy for dflow-extender in dflow is None, we actually change the image-pull policy to "IfNotPresent" in rid-kit, so that once you load the image to the minikube, the pods will directly use it without pulling it.

**You may encouter this message during the workflow**
```
unschedulable: 0/1 nodes are available: 1 too many pods. preemption: 0/1 nodes are available: 1 no preemption victims found for incoming pod.
```
Usually this is not a big deal, the reason for this warning is that the running pods exceed the maximum number of pods setting in k8s. The default setting for max number of pods is 110, this is a reasonable number. If you really want to enlarge the max number of pods, there are two ways of doing it, you can set it when starting the minikube by
```
minikube start --extra-config=kubelet.max-pods=200
```
this set the max pods number to 200. If you do not want to restart the minikube, you can also change the setting by entering the minikube node, and manually change the settings in var/lib/kubelet/config.yaml (If it is not the right path type "ps -ef | grep -i kubelet" to check it).

**In rare cases you may encouter this problem after submitting the workflow**
```
NAME                        STATUS                AGE   DURATION   PRIORITY
reinforced-dynamics-xkghb   Pending               30m   0s        0
```
namely the workflow gets stuck in pending status for a long time. This problem is probably due to an update of the argo image, you 
may check "argoproj/workflow-controller" or "argoproj/argocli" image on dockerhub to confirm. When this is actually the problem, you can use the stable version of yaml file in dflow repo, which is located at "dflow/manifests/quick-start-postgres-stable-cn.yaml".