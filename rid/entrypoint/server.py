import os, sys, logging
from rid.utils import run_command
from rid.constants import argo_namespace

logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

PORTS = {
    "argo": 2746,
    "minio-server": 9000,
    "minio-ui": 9001
}


def check_port_status():
    return_code, out, err = run_command(['ps', 'aux'])
    assert return_code == 0
    job_list = out.split("\n")
    port_status = {
        "argo": 0,
        "minio-server": 0,
        "minio-ui": 0,
    }
    for job in job_list:
        if "port-forward" in job and "kubectl" in job:
            for server in port_status.keys():
                port = PORTS[server]
                if f"{port}:{port}" in job:
                    port_status[server] = 1
    return port_status


def forward_ports():
    port_status = check_port_status()
    if port_status["argo"] == 1:
        logger.info('Port "agro-server" has been launched and running.')
    else:
        port = PORTS["argo"]
        os.system(f"nohup minikube kubectl -- -n {argo_namespace} port-forward deployment/argo-server {port}:{port} --address 0.0.0.0 > /tmp/argo-server.out 2>&1 &")
    
    if port_status["minio-server"] == 1:
        logger.info('Port "minio-server" has been launched and running.')
    else:
        port = PORTS["minio-server"]
        os.system(f"nohup minikube kubectl -- -n {argo_namespace} port-forward deployment/minio {port}:{port} --address 0.0.0.0 > /tmp/minio-server.out 2>&1 &")
    
    if port_status["minio-ui"] == 1:
        logger.info('Port "minio-ui" has been launched and running.')
    else:
        port = PORTS["minio-ui"]
        os.system(f"nohup minikube kubectl -- -n {argo_namespace} port-forward deployment/minio {port}:{port} --address 0.0.0.0 > /tmp/minio-ui.out 2>&1 &")
