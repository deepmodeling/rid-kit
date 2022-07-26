import argparse, os, json, glob, sys, logging
from pathlib import Path

from dflow import (
    Workflow,
    Step,
    Steps,
    upload_artifact,
    download_artifact,
)
from typing import (
    Optional,
    List,
)
from .submit import submit_rid
from .info import information
from .server import forward_ports
from rid import (
    __version__
)

logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def main_parser() -> argparse.ArgumentParser:
    """RiD commandline options argument parser.
    Notes
    -----
    This function is used by documentation.
    Returns
    -------
    argparse.ArgumentParser
        the argument parser
    """
    parser = argparse.ArgumentParser(
        # description="RiD: Enhanced sampling methods under the concurrent learning framework.",
        description=information,
        # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n"
    )
    subparsers = parser.add_subparsers(title="Valid subcommands", dest="command")

    parser_run = subparsers.add_parser(
        "port-forward",
        help="Port forward",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="\n"
    )
    
    parser_run = subparsers.add_parser(
        "status",
        help="Task status and Server status",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="\n"
    )
    
    # workflow submit.
    parser_run = subparsers.add_parser(
        "submit",
        help="Submit RiD workflow",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="\n"
    )
    parser_run.add_argument(
        "--mol", "-i", help="Initial conformation files.", dest="mol",
    )
    # parser_run.add_argument(
    #     "--conf", "-i", help="Initial conformation files.", dest="model"
    # )
    # parser_run.add_argument(
    #     "--top", "-p", help="Topology file.", dest="top"
    # )
    # parser_run.add_argument(
    #     "--forcefield", "-f", help="Forcefield files.", dest="ff"
    # )
    parser_run.add_argument(
        "--config", "-c", help="RiD configuration.", dest="config"
    )
    parser_run.add_argument(
        "--machine", "-m", help="Machine configuration.", dest="machine"
    )
    
    # Explore
    parser_exp = subparsers.add_parser(
        "mdrun",
        help="Submit RiD workflow",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_exp.add_argument(
        "--mol", "-i", help="Initial conformation files.", dest="mol",
    )
    # parser_exp.add_argument(
    #     "--conf", "-i", help="Initial conformation files.", dest="model"
    # )
    # parser_exp.add_argument(
    #     "--top", "-p", help="Topology file.", dest="top"
    # )
    # parser_exp.add_argument(
    #     "--forcefield", "-f", help="Forcefield files.", dest="ff"
    # )
    parser_exp.add_argument(
        "--config", "-c", help="RiD configuration.", dest="config"
    )
    parser_exp.add_argument(
        "--machine", "-m", help="Machine configuration.", dest="machine"
    )
    
    # train
    parser_train = subparsers.add_parser(
        "train",
        help="Train RiD neural networks",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_train.add_argument(
        "--data", "-d", help="Training data."
    )
    parser_train.add_argument(
        "--config", "-c", help="RiD configuration."
    )

    # NN dimension reduction.
    parser_redim = subparsers.add_parser(
        "redim",
        help="NN dimension reduction by Monte Carlo integral.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_redim.add_argument(
        "--networks", "-n", help="Training data."
    )
    parser_redim.add_argument(
        "--dim1", help="Dimension 1."
    )
    parser_redim.add_argument(
        "--dim2", help="Dimension 2."
    )

    # --version
    parser.add_argument(
        '--version', 
        action='version', 
        version='RiD v%s' % __version__,
    )

    return parser


def parse_args(args: Optional[List[str]] = None):
    """RiD commandline options argument parsing.
    Parameters
    ----------
    args: List[str]
        list of command line arguments, main purpose is testing default option None
        takes arguments from sys.argv
    """
    parser = main_parser()

    parsed_args = parser.parse_args(args=args)
    if parsed_args.command is None:
        parser.print_help()

    return parsed_args
    

def main():
    args = parse_args()
    logger.info("{}\n{}".format("Software Statement", information))
    # print(args)
    if args.command == "submit":
        
        logger.info("Preparing RiD ...")
        mol_path = Path(args.mol)
        confs = glob.glob(str(mol_path.joinpath("*.gro")))
        assert len(confs) > 0, "No valid conformation files."
        top = glob.glob(str(mol_path.joinpath("*.top")))
        assert len(top) > 0, "No valid topology files."
        top_file = top[0]
        forcefield = glob.glob(str(mol_path.joinpath("*.ff")))
        if len(forcefield) == 0:
            forcefield = None
        else:
            forcefield = forcefield[0]
        models = glob.glob(str(mol_path.joinpath("*.pb")))
        if len(models) == 0:
            models = None
        submit_rid(
            confs = confs,
            topology = top_file,
            rid_config = args.config,
            machine_config = args.machine,
            models = models,
            forcefield = forcefield
        )
        logger.info('The task is displayed on "https://127.0.0.1:2746".')
        logger.info('Artifacts (Files) are listed on "https://127.0.0.1:9001".')
    elif args.command == "explore":
        logger.info("RiD Exploration.")
        return None
    elif args.command == "port-forward":
        forward_ports()
    else:
        raise RuntimeError(f"unknown command {args.command}")