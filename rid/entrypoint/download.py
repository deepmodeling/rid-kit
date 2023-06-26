from dflow import Workflow, download_artifact
import numpy as np
import os
from datetime import datetime

from dflow import config, s3_config
from dflow.plugins import bohrium
from dflow.plugins.bohrium import TiefblueClient


def download_rid(workflow_id, pod, file, iteration_start, iteration_end, outputs):
    wf = Workflow(id=workflow_id)
    all_steps = wf.query_step()

    iter_a = int(iteration_start)
    iter_e = int(iteration_end)+1
    print("start from iteration",iter_a)
    print("end to iteration",iter_e-1)
    for step in all_steps:
        if step["type"] == "Pod":
            for index in range(iter_a, iter_e):
                # print(step["key"])
                if step['key'] == "iter-%03d-%s-0"%(index,pod):
                    print("downloading "+step["key"])
                    download_artifact(step.outputs.artifacts[file], path="%s/iter%s"%(outputs,index))
                elif step['key'] == "iter-%03d-%s"%(index,pod):
                    print("downloading "+step["key"])
                    download_artifact(step.outputs.artifacts[file], path="%s/iter%s"%(outputs,index))
                elif step["key"] == "000-%s-0"%pod:
                    print("downloading "+step["key"])
                    download_artifact(step.outputs.artifacts[file], path="%s"%outputs)
                    break
                elif step["key"] == "000-%s"%pod:
                    print("downloading "+step["key"])
                    download_artifact(step.outputs.artifacts[file], path="%s"%outputs)
                    break

# with open("time_profile.txt","w") as f:
#     for step in all_steps:
#         if step['key'] is not None:
#             f.write("step key:"+str(step["key"]))
#         if step['startedAt'] is not None:
#             start_time_str = str(step['startedAt'])
#             f.write(" start time:"+str(start_time_str))
#             start_time = datetime.strptime(start_time_str, "%Y-%m-%dT%H:%M:%SZ")
#         if step['finishedAt'] is not None:
#             end_time_str = str(step['finishedAt'])
#             f.write(" finish time:"+str(end_time_str))
#             end_time = datetime.strptime(end_time_str, "%Y-%m-%dT%H:%M:%SZ")
#             # Calculate the duration of the task
#             duration = end_time - start_time
#             f.write(" duration:"+str(duration))
#         f.write("\n")