import argparse, os, glob, sys, logging
from pathlib import Path
from typing import (
    Optional,
    List,
)


NUMEXPR_MAX_THREADS = os.getenv("NUMEXPR_MAX_THREADS")
if NUMEXPR_MAX_THREADS is None:
    NUMEXPR_MAX_THREADS = 8
    os.environ["NUMEXPR_MAX_THREADS"] = str(NUMEXPR_MAX_THREADS)

try:
    import tensorflow.compat.v1 as tf
    tf.logging.set_verbosity(tf.logging.ERROR)
    tf.disable_v2_behavior()
except ImportError:
    import tensorflow as tf
    tf.logging.set_verbosity(tf.logging.ERROR)


from .submit import submit_rid
from .info import information
from .server import forward_ports
from .cli import rid_ls, rid_rm
from rid import __version__


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
    
    subparsers.add_parser(
        "ls",
        help="List all rid tasks.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers.add_parser(
        "dp",
        help="Something interesting.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser_rm = subparsers.add_parser(
        "rm",
        help="Remove workflows of RiD.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_rm.add_argument(
        "WORKFLOW_ID", help="Workflow ID"
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
    parser_run.add_argument(
        "--config", "-c", help="RiD configuration.", dest="config"
    )
    parser_run.add_argument(
        "--machine", "-m", help="Machine configuration.", dest="machine"
    )

    # resubmit
    parser_rerun = subparsers.add_parser(
        "resubmit",
        help="Submit RiD workflow",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_rerun.add_argument(
        "ID", help="Workflow ID."
    )
    parser_rerun.add_argument(
        "--mol", "-i", help="Initial conformation files.", dest="mol",
    )
    parser_rerun.add_argument(
        "--config", "-c", help="RiD configuration.", dest="config"
    )
    parser_rerun.add_argument(
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
    elif args.command == "ls":
        rid_ls()
    elif args.command == "rm":
        rid_rm(args.WORKFLOW_ID)
    elif args.command == "dp":
        logger.info("Molecule Simulates the Future!")
    else:
        raise RuntimeError(f"unknown command {args.command}")