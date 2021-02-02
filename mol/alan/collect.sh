#!/bin/bash

if test $# -ne 1; then
    echo usage 
    echo $0 src 
    exit
fi

src=$1
if test ! -d $src; then
    echo no dir $src
    exit
fi

cwd=`pwd`
cd $src
src=`pwd`
cd $cwd

ffield=`grep pdb2gmx $src/01.pdb2gmx/gmx_pdb2gmx.log | awk '{for(i=1;i<=NF;i++)if($i~/-ff/)print $(i+1)}'`
nalpha=`grep CA $src/mol.pdb | wc -l`
pa=`printf %02d $nalpha`

echo nalpha: $pa
echo force field: $ffield

if test -d $ffield/$pa.gas || test -d $ffield/$pa.sol; then
    echo existing folder $ffield/$pa.gas or $ffield/$pa.sol, do nothing
    exit
fi

if test ! -d $ffield; then
    mkdir $ffield
    cd $ffield 
    ln -s ../grompp.gas.mdp .
    ln -s ../grompp.mdp .
    cd ..
fi

mkdir $ffield/$pa.gas
cd $ffield/$pa.gas
cp -L $src/03.regzt/confout.gro ./conf.gro
cp -L $src/03.regzt/topol.top ./topol.top
ln -s ../grompp.gas.mdp ./grompp.mdp
cd ../..

mkdir $ffield/$pa.sol
cd $ffield/$pa.sol
cp -L $src/05.equi/confout.gro ./conf.gro
cp -L $src/05.equi/topol.top ./topol.top
ln -s ../grompp.mdp ./grompp.mdp
cd ../
ln -sf $pa.sol $pa
cd ..
