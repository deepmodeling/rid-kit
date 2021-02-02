#!/bin/bash

./tools/pbc.angle.py plm.res.out pbc.res.out

nline=`grep -v \# pbc.res.out | wc -l`
nword=`grep -v \# pbc.res.out | wc -w`
ncloumn=`echo "$nword / $nline " | bc`
clist=""
for ii in `seq 2 $ncloumn`
do
    if test $ii -eq 2; then
	clist=$ii
    else
	clist="$clist,$ii"
    fi
done
# echo $clist

msp_avg -f pbc.res.out -m $clist -n 8 -t .9 | grep -v \# > avgins.centers.out

grep KAPPA plumed.res.dat  | awk '{print $4}' | cut -d '=' -f 2 > kappa.out

./tools/cmpf.py
# python -c "import numpy as np; kk=np.loadtxt('kappa.out'); cc=np.loadtxt('centers.out'); data=np.loadtxt ('avgins.centers.out'); avgins=data[:,0]; ff=np.multiply (-kk , (avgins - cc) ); np.savetxt ('force.out', ff, fmt='%.10e'); erravg = data[:,1]; fe = np.multiply (kk, erravg); np.savetxt ('ferror.out', fe, fmt='%.10e');"

cat force.out | tr '\n' ' ' > tmp.out
mv -f tmp.out force.out
echo "" >> force.out

cat ferror.out | tr '\n' ' ' > tmp.out
mv -f tmp.out ferror.out
echo "" >> ferror.out

