import argparse, os, glob, sys, logging
from operator import index
from pathlib import Path
from typing import (
    Optional,
    List,
)
import rid

if os.getenv('DFLOW_DEBUG'):
    from dflow import config
    config["mode"] = "debug"
    
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
from .resubmit import resubmit_rid
from .explore import explore_rid
from .label import label_rid
from .relabel import relabel_rid
from .redim import redim_rid
from .reredim import reredim_rid
from .train import train_rid
from .download import download_rid
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
    parser_run.add_argument(
        "--id", "-d", help="Workflow user defined ID", default = None, dest="id"
    )

    # resubmit
    parser_rerun = subparsers.add_parser(
        "resubmit",
        help="Resubmit RiD workflow",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_rerun.add_argument(
        "WORKFLOW_ID", help="Workflow ID."
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
    parser_rerun.add_argument(
        "--id", "-d", help="Workflow user defined ID", default = None, dest="id"
    )
    parser_rerun.add_argument(
        "--iteration", "-t", help="restart from t-th iteration.", default = None, dest="iteration"
    )
    parser_rerun.add_argument(
        "--pod", "-p", help="restart from the pod.", default = None, dest="pod"
    )
    
    # Explore
    parser_exp = subparsers.add_parser(
        "explore",
        help="Exploration MD",
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
    
    # Label
    parser_label = subparsers.add_parser(
        "label",
        help="labeling MD workflow",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_label.add_argument(
        "--mol", "-i", help="Initial conformation files.", dest="mol",
    )
    parser_label.add_argument(
        "--config", "-c", help="RiD configuration.", dest="config"
    )
    parser_label.add_argument(
        "--machine", "-m", help="Machine configuration.", dest="machine"
    )
    
    # Relabel
    parser_relabel = subparsers.add_parser(
        "relabel",
        help="Resubmit labeling MD workflow",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_relabel.add_argument(
        "WORKFLOW_ID", help="Workflow ID."
    )
    parser_relabel.add_argument(
        "--mol", "-i", help="Initial conformation files.", dest="mol",
    )
    parser_relabel.add_argument(
        "--config", "-c", help="RiD configuration.", dest="config"
    )
    parser_relabel.add_argument(
        "--machine", "-m", help="Machine configuration.", dest="machine"
    )
    
    # train
    parser_train = subparsers.add_parser(
        "train",
        help="Train RiD neural networks",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_train.add_argument(
        "--mol", "-i", help="Training data.", dest="mol",
    )
    parser_train.add_argument(
       "--config", "-c", help="RiD configuration.", dest="config"
    )
    parser_train.add_argument(
         "--machine", "-m", help="Machine configuration.", dest="machine"
    )

    # NN dimension reduction.
    parser_redim = subparsers.add_parser(
        "redim",
        help="NN dimension reduction by Monte Carlo integral.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_redim.add_argument(
        "--mol", "-i", help="Neural networks path", dest="mol",
    )
    parser_redim.add_argument(
        "--config", "-c", help="RiD configuration.", dest="config"
    )
    parser_redim.add_argument(
        "--machine", "-m", help="Machine configuration.", dest="machine"
    )
    parser_redim.add_argument(
        "--id", "-d", help="Workflow user defined ID", default = None, dest="id"
    )
    
    # resubmit NN dimension reduction.
    parser_reredim = subparsers.add_parser(
        "reredim",
        help="NN dimension reduction by Monte Carlo integral.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_reredim.add_argument(
        "--mol", "-i", help="Neural networks path", dest="mol",
    )
    parser_reredim.add_argument(
        "--config", "-c", help="RiD configuration.", dest="config"
    )
    parser_reredim.add_argument(
        "--machine", "-m", help="Machine configuration.", dest="machine"
    )
    parser_reredim.add_argument(
        "--id", "-d", help="Workflow user defined ID", default = None, dest="id"
    )
    parser_reredim.add_argument(
        "WORKFLOW_ID", help="Workflow ID."
    )
    parser_reredim.add_argument(
        "--pod", "-p", help="restart from the pod.", default = None, dest="pod"
    )
    # download
    parser_download = subparsers.add_parser(
        "download",
        help="Download files in RiD workflow",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_download.add_argument(
        "WORKFLOW_ID", help="Workflow ID."
    )
    parser_download.add_argument(
        "--pod", "-p", help="The pod in RiD workflow", dest="pod"
    )
    parser_download.add_argument(
        "--file", "-f", help="The file to download inside the step", dest = "file"
    )
    parser_download.add_argument(
        "--iteration_start", "-a", help="download from a-th iteration.", default = 1, dest="iteration_start"
    )
    parser_download.add_argument(
        "--iteration_end", "-e", help="download to e-th iteration.", default = 100, dest="iteration_end"
    )
    parser_download.add_argument(
        "--outputs", "-o", help="output directory of downloading", default = "./", dest="outputs"
    )

    # --version
    parser.add_argument(
        '--version', 
        action='version', 
        version='RiD v%s' % __version__,
    )

    return parser


def parse_args(args: Optional[List[str]] = None):
    """
    RiD commandline options argument parsing.
    
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


def parse_submit(args):
    mol_path = Path(args.mol)
    allfiles = glob.glob(str(mol_path.joinpath("*")))
    confs = []
    top_file = None
    forcefield = None
    index_file = None
    plm_files = []
    dp_files = []
    models = []
    data_file = None
    otherfiles = []
    for file in allfiles:
        if os.path.basename(file).endswith("ff"):
            forcefield = file
        elif os.path.basename(file).endswith("gro") or os.path.basename(file).endswith("lmp"):
            confs.append(file)
        elif os.path.basename(file).endswith("top"):
            top_file = file
        elif os.path.basename(file).endswith("pb"):
            models.append(file)
        elif os.path.basename(file).endswith("npy"):
            data_file = file
        elif os.path.basename(file).endswith("ndx"):
            index_file = file
        elif os.path.basename(file).endswith("json") or os.path.basename(file).endswith("raw"):
            dp_files.append(file)
        elif os.path.basename(file).endswith("out"):
            plm_files.append(file)
        else:
            otherfiles.append(file)
        
    return confs, top_file, models, forcefield, index_file, data_file, dp_files, plm_files, otherfiles


def log_ui():
    if os.getenv("DFLOW_HOST") is not None:
        logger.info('The task is displayed on %s.'%os.getenv("DFLOW_HOST"))
    else:
        logger.info('The task is displayed on "https://127.0.0.1:2746".')
        logger.info('Artifacts (Files) are listed on "https://127.0.0.1:9001".')


def main():
    args = parse_args()
    if args.command == "submit":
        logger.info("Preparing RiD ...")
        confs, top_file, models, forcefield, index_file, data_file, dp_files,plm_files, otherfiles = parse_submit(args)
        submit_rid(
            confs = confs,
            topology = top_file,
            rid_config = args.config,
            machine_config = args.machine,
            workflow_id_defined = args.id,
            models = models,
            forcefield = forcefield,
            index_file = index_file,
            data_file = data_file,
            dp_files = dp_files,
            otherfiles = otherfiles
        )
        log_ui()
    elif args.command == "resubmit":
        logger.info("Preparing RiD ...")
        confs, top_file, models, forcefield, index_file, data_file, dp_files,plm_files, otherfiles = parse_submit(args)
        resubmit_rid(
            workflow_id=args.WORKFLOW_ID,
            confs = confs,
            topology = top_file,
            rid_config = args.config,
            machine_config = args.machine,
            workflow_id_defined = args.id,
            iteration = args.iteration,
            pod = args.pod,
            models = models,
            forcefield = forcefield,
            index_file = index_file,
            data_file = data_file,
            dp_files = dp_files,
            otherfiles = otherfiles
        )
        log_ui()
    elif args.command == "explore":
        logger.info("RiD Exploration ..")
        confs, top_file, models, forcefield, index_file, data_file, dp_files,plm_files, otherfiles = parse_submit(args)
        explore_rid(
            confs = confs,
            topology = top_file,
            rid_config = args.config,
            machine_config = args.machine,
            models = models,
            forcefield = forcefield,
            index_file = index_file,
            dp_files = dp_files,
            otherfiles = otherfiles
        )
        log_ui()
    elif args.command == "label":
        logger.info("Labeling MD ...")
        confs, top_file, models, forcefield, index_file, data_file, dp_files,plm_files, otherfiles = parse_submit(args)
        label_rid(
            confs = confs,
            topology = top_file,
            rid_config = args.config,
            machine_config = args.machine,
            models = models,
            forcefield = forcefield,
            index_file = index_file,
            dp_files = dp_files,
            otherfiles = otherfiles
        )
        log_ui()
    elif args.command == "relabel":
        logger.info("Labeling MD ...")
        confs, top_file, models, forcefield, index_file, data_file, dp_files,plm_files, otherfiles = parse_submit(args)
        relabel_rid(
            workflow_id=args.WORKFLOW_ID,
            confs = confs,
            topology = top_file,
            rid_config = args.config,
            machine_config = args.machine,
            models = models,
            forcefield = forcefield,
            index_file = index_file,
            dp_files = dp_files,
            otherfiles = otherfiles
        )
        log_ui()
    elif args.command == "redim":
        logger.info("Preparing MCMC ...")
        confs, top_file, models, forcefield, index_file, data_file, dp_files,plm_files, otherfiles = parse_submit(args)
        redim_rid(
            rid_config = args.config,
            machine_config = args.machine,
            models = models,
            plm_out = plm_files,
            workflow_id_defined = args.id,
        )
        log_ui()
    elif args.command == "reredim":
        logger.info("Preparing MCMC ...")
        confs, top_file, models, forcefield, index_file, data_file, dp_files,plm_files, otherfiles = parse_submit(args)
        reredim_rid(
            workflow_id=args.WORKFLOW_ID,
            rid_config = args.config,
            machine_config = args.machine,
            workflow_id_defined = args.id,
            models = models,
            plm_out = plm_files,
            pod = args.pod
        )
        log_ui()
    elif args.command == "train":
        logger.info("Training RiD FES ...")
        confs, top_file, models, forcefield, index_file, data_file, dp_files,plm_files, otherfiles = parse_submit(args)
        train_rid(
            data = data_file,
            rid_config = args.config,
            machine_config = args.machine
        )
    elif args.command == "download":
        logger.info("Downloading files ...")
        download_rid(
            workflow_id = args.WORKFLOW_ID,
            pod = args.pod,
            file = args.file,
            iteration_start = args.iteration_start,
            iteration_end = args.iteration_end,
            outputs = args.outputs
        )
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