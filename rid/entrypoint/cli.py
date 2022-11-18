import logging, sys, os
from rid.utils import run_command
from rid.constants import argo_namespace


logging.basicConfig(
    format="%(asctime)s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def rid_ls():
    return_code, out, err = run_command(["argo", "list", "-n", argo_namespace])
    assert return_code == 0, err
    logger.info(f"\n\n\tReinforced Dynamics Workflow\n\n{out}")

def rid_rm(workflow_id):
    return_code, out, err = run_command(["argo", "delete", "-n", argo_namespace, workflow_id])
    assert return_code == 0, err
    logger.info(f"\n{out}")

