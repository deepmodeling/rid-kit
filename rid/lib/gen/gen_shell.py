import os

def gen_res_shell(shell_path):
    ret = """#!/bin/bash

targets=$*
res_names=`grep RESTRAINT plumed.res.templ | cut -d ':' -f 1`
nres=`echo $res_names | wc -w`

if test $nres -ne $#; then
    echo input number of angles $# does not match number of restraint $nres
    exit 
fi    

rm -f plumed.res.dat
cp plumed.res.templ plumed.res.dat

count=0
centers=""

for ii in $res_names
do
    arad=$1
    centers="$centers $arad"
    sed -i "/$ii/s/AT=.*/AT=$arad/g" plumed.res.dat
    shift 1
done

echo $centers > centers.out
"""
    if os.path.basename(shell_path) == '':
        shell_path = os.path.abspath(shell_path) + "/general_mkres.sh"
    with open(shell_path, 'w') as sh:
        sh.write(ret)
    print(shell_path + " was generated.")

if __name__ == '__main__':
    gen_res_shell("./")
