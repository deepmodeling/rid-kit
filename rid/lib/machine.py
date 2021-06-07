import json

from dpdispatcher.lazy_local_context import LazyLocalContext
from dpdispatcher.submission import Resources
from dpdispatcher.pbs import PBS
from dpdispatcher.slurm import Slurm
from dpdispatcher.machine import Machine


def set_resource(json_file, target='enhcMD'):
    fp = open(json_file, 'r')
    jdata = json.load(fp)
    fp.close()
    resources_dict = jdata[target]["resources"]
    resources = Resources.load_from_dict(resources_dict)
    return resources


def set_machine(json_file, target="enhcMD"):
    fp = open(json_file, 'r')
    jdata = json.load(fp)
    fp.close()
    machine_dict = jdata[target]["machine"]
    machine = Machine.load_from_dict(machine_dict)
    # machine_dict["batch_type"] = batch_dict["batch_type"].lower()
    return machine


if __name__ == "__main__":
    set_machine("/home/dongdong/wyz/rid-kit/examples/machine.json", "cmpf")
