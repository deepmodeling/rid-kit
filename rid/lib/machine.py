import json

from dpdispatcher.lazy_local_context import LazyLocalContext
from dpdispatcher.submission import Resources
from dpdispatcher.pbs import PBS
from dpdispatcher.slurm import Slurm
from dpdispatcher.batch_object import BatchObject, Machine

def set_resource(json_file, target='enhcMD'):
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    fp.close()
    resource_dict = jdata[target]["resources"]
    return Resources(**resource_dict) 

def set_batch(json_file, target="enhcMD"):
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    fp.close()
    batch_dict = jdata[target]["batch"]
    batch_dict["batch_type"] = batch_dict["batch_type"].lower()
    return BatchObject(jdata=batch_dict)

    

if __name__ == "__main__":
    set_batch("/home/dongdong/wyz/rid-kit/examples/machine.json", "cmpf")