import subprocess
from typing import Optional

def run_command(
        cmd: str, 
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
        if stdin is not None:
            subp.stdin.write(stdin)
            subp.stdin.close()
        subp.wait()
        out, err = subp.communicate()
        return_code = subp.poll()
    return return_code, out, err

