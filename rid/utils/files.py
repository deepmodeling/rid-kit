import json
from typing import Dict


def write_txt(
        fname: str, 
        fcont: str
    ):
    with open(fname, "w") as fn:
        fn.write(fcont)


def load_json(
        fname: str
    ) -> Dict:
    with open(fname, "r") as fn:
        jdata = json.load(fn)
    return jdata


def dump_json(
        fname: str, 
        fcont: Dict,
        indent: int = 4
    ):
    with open(fname, "w") as fn:
        json.dump(fcont, fn, indent=indent)