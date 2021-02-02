#!/bin/bash

targets=$*
# kk=500
nalpha=10
nangle=4
alpha_idx_fmt=%02d
angle_idx_fmt=%02d

# echo num args: $#
tot_n=`echo "$nangle * $nalpha" | bc `
if test $# -ne $tot_n; then
    echo should input $tot_nangle angle in the unit of rad
    exit
fi

rm -f plumed.res.dat
cp plumed.res.templ plumed.res.dat

count=0
centers=""
for ii in $targets; 
do
#    arad=`echo "$ii / 180. * 3.141593" | bc -l`
    arad=$ii
    centers="$centers $arad"
    ialpha=`echo "$count / $nangle" | bc`
    pial=`printf $alpha_idx_fmt $ialpha`
    iangle=`echo "$count % $nangle" | bc `    
    pian=`printf $angle_idx_fmt $iangle`
    
    sed -i "/res-$pial-$pian/s/AT=.*/AT=$arad/g" plumed.res.dat

    count=$(($count+1))
done

# sed -i "s/KAPPA=.* /KAPPA=$kk /g" plumed.res.dat

echo $centers > centers.out
