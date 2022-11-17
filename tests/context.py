import sys,os
rid_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, rid_path)
import rid

if os.getenv('SKIP_UT_WITH_DFLOW'):
    skip_ut_with_dflow = (int(os.getenv('SKIP_UT_WITH_DFLOW')) != 0)
    skip_ut_with_dflow_reason = 'skip because environment variable SKIP_UT_WITH_DFLOW is set to non-zero'
else:
    skip_ut_with_dflow = False
    skip_ut_with_dflow_reason = ''

# one needs to set proper values for the following variable.
default_image = 'pkufjhdocker/rid-tf-cpu:latest'
default_host = None
if os.getenv('DFLOW_DEBUG'):
    from dflow import config
    config["mode"] = "debug"