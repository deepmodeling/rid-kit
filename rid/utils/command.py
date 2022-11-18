import subprocess
from typing import Optional, List

def run_command(
        cmd: List, 
        stdin: Optional[str] = None,
        shell: Optional[bool] = None
    ):

    stdout = subprocess.PIPE
    stderr = subprocess.PIPE
    
    with subprocess.Popen(
        args=cmd,
        stdin=subprocess.PIPE,
        stdout=stdout,
        stderr=stderr,
        encoding="utf-8",
        shell=shell
    ) as subp:
        out, err = subp.communicate(input=stdin)
        return_code = subp.poll()
    return return_code, out, err
