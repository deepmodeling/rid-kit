import os, glob, json
from ridkit.lib.utils import create_path

def gen_rid (out_dir, mol_dir, rid_json) :
    mol_dir = os.path.abspath(mol_dir) + "/"
    out_dir = os.path.abspath(out_dir) + "/"
    conf_file = glob.glob(mol_dir + "*.gro")
    rid_json = os.path.abspath(rid_json)
    fp = open (rid_json, 'r')
    jdata = json.load (fp)
    numb_walkers = jdata['numb_walkers']
    create_path(out_dir)
    assert (len(conf_file) >= int(numb_walkers)), "number of conformation files must be equal to the number of walkers."
    print("RiD dir has prepared.")


def main_debug():
    create_path("./debug")
    return
   
if __name__ == '__main__':
    main_debug()
