import os, sys
import logging
from typing import Sequence
from rid.common.gromacs.gmx_constant import gmx_trjconv_cmd, gmx_traj_cmd
from rid.constants import gmx_coord_name, gmx_force_name
from rid.utils import list_to_string
from rid.utils import run_command


logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

def generate_coords(
    system: str,
    trr: str,
    top: str,
    out_coord: str = gmx_coord_name
):
    echo_string = "echo -e '%s\n' | "%system
    cmd_list = gmx_traj_cmd.split(" ")
    cmd_list += ["-f", str(trr)]
    cmd_list += ["-s", str(top)]
    cmd_list += ["-ox", out_coord]
    os.system(echo_string+list_to_string(cmd_list, " "))
    
def generate_forces(
    system: str,
    trr: str,
    top: str,
    out_force: str = gmx_force_name
):
    echo_string = "echo -e '%s\n' | "%system
    cmd_list = gmx_traj_cmd.split(" ")
    cmd_list += ["-f", str(trr)]
    cmd_list += ["-s", str(top)]
    cmd_list += ["-of", out_force]
    os.system(echo_string+list_to_string(cmd_list, " "))

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
    

    