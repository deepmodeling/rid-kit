import os, sys
import logging
from typing import Sequence
from rid.common.gromacs.gmx_constant import gmx_trjconv_cmd
from rid.utils import list_to_string
from rid.utils import run_command


logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def slice_trjconv(
        xtc: str,
        top: str,
        selected_time: float,
        output_group: int = 0,
        output: str = "conf.gro"
    ):
    logger.info("slicing trajectories by gmx trjconv command ...")
    logger.warning("You are using `gmx trjconv` to slice trajectory, "
    "make sure that the selected index is in the unit of time (ps).")
    cmd_list = gmx_trjconv_cmd.split()
    cmd_list += ["-f", str(xtc)]
    cmd_list += ["-s", str(top)]
    cmd_list += ["-dump", str(selected_time)]
    cmd_list += ["-o", output]
    logger.info(list_to_string(cmd_list, " "))
    return_code, out, err = run_command(
        cmd_list,
        stdin=f"{output_group}\n"
    )
    assert return_code == 0, err
    

    