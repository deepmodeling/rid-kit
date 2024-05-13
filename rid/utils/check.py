# Trouble shooting for common errors during the workflow

import os
import subprocess
import numpy as np
from pathlib import Path

def check_cv_file(file_list : list):
    """Parse the cv files with plumed to check if the files are valid.
    
    Parameters
    ----------
    file_list : list
        List of file absolute paths. Only one .pdb file is allowed.
        
    Returns
    -------
    Rsl : bool
        True if the file is valid, False otherwise.
    Rmsg : str
        The whole parsing message of PLUMED if the file is invalid, 
        none otherwise.
    """
    input_dir = Path(file_list[0]).parent
    os.chdir(input_dir)
    
    cv_file_list = [file for file in file_list if file[-4:] != ".pdb"]
    strut_pdb = [file for file in file_list if file[-4:] == ".pdb"]
    assert len(strut_pdb) == 1, \
        "There should be only one .pdb file in the cv files."
    
    # Fetch the number of atoms from the .pdb file
    with open(strut_pdb[0], "rb") as f:
        block = -1024
        flag = True
        while flag:
            f.seek(block, 2)
            lines = f.readlines()
            lines.reverse()
            for line in lines:
                if line.startswith(b'ATOM'):
                    natoms = int(line[6:11])
                    flag = False
                    break
            block *= 2
    
    # Parse the file with plumed
    Rsl = True
    Rmsg = ''
    for file in cv_file_list:
        cmd = f"plumed driver --plumed {file} \
            --parse-only --natoms {natoms}"
        (Status, Rmsg) = subprocess.getstatusoutput(cmd)
        if Status != 0:
            Rsl = False
            break
        
    return Rsl, Rmsg
