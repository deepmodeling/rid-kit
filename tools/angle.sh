#!/bin/bash

set -eu

nalpha=10
nangle=4

rm -f angdist.*.xvg
rm -f angaver.*.xvg

angcmd="gmx angle"
if test -f traj_comp.xtc; then
    file=traj_comp.xtc
else
    file=traj.xtc
fi
if test -f traj.trr; then
    file=traj.trr
fi

nalpha_=`echo "$nalpha - 1" | bc`
nalpha_=`printf %02d $nalpha_`
nangle_=`echo "$nangle - 1" | bc`
nangle_=`printf %02d $nangle_`

for ii in `seq -w 00 $nalpha_`; 
do
    for jj in `seq -w 00 $nangle_`;
    do
	code=`echo "$ii * $nangle + $jj" | bc `	
	echo doing $ii $jj
	echo $code | $angcmd -f $file -n angle.ndx  -type dihedral -ov angaver.$ii.$jj.xvg -od angdist.$ii.$jj.xvg -xvg none &> /dev/null
	awk '{print $2}' angaver.$ii.$jj.xvg > tmp.$ii.$jj.xvg
    done
    paste tmp.$ii.*.xvg > angle.deg.$ii.out
    rm -f tmp.$ii.*.xvg
    python3 -c "import numpy as np; d = np.loadtxt (\"angle.deg.$ii.out\"); np.savetxt (\"angle.rad.$ii.out\", np.reshape(d *0.017453292519943295, [-1,$nangle]), fmt=\"%.6f\")"
done

paste angle.rad.*.out | tr '\t' ' ' > angle.rad.out
paste angle.deg.*.out | tr '\t' ' ' > angle.deg.out
