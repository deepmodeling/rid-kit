#!/bin/bash

set -e

if test $# -ne 3; then
    echo usage
    echo $0 forcefield nalpha box
    exit
fi

echo using ff $1 box $3 nalpha $2
forcefield=$1
box=$3
nalpha=$2
cwd=`pwd`
plm_tool_dir=$HOME/study/deep.fe/source/ala-n/tools
gmx_dir=$HOME/local/gromacs/5.1.4-plumed-fe/
source $gmx_dir/bin/GMXRC

function cut_pdb_11 () {
    nalpha=$nalpha
    file=alphaHelix.pdb
    cp $file mol.pdb
    a_start=`echo "$nalpha + 2" | bc `
    a_end=11
    for ii in `seq $a_start $a_end`;
    do
	grep -v "ALA  *$ii" mol.pdb > tmp.out
	mv -f tmp.out mol.pdb
    done    
    if test $a_start -ge 10; then
	sed -i "/NME/s/12/$a_start/g" mol.pdb
    else
	sed -i "/NME/s/12/ $a_start/g" mol.pdb
    fi
}

function f01_pdb2gmx () {
    dir=01.pdb2gmx
    if test -d $dir; then
	echo existing dir $dir, do nothing.
	exit
    fi
    mkdir $dir
    cd $dir/
    ln -s ../mol.pdb .
    gmx pdb2gmx -f mol.pdb -ff $forcefield -water tip3p &> gmx_pdb2gmx.log

    natoms=`cat conf.gro | head -n 2 | tail -n 1`
    nhead=`echo "$natoms + 2" | bc`
    head -n $nhead conf.gro > tmp.gro
    echo "$box $box $box" >> tmp.gro
    mv -f tmp.gro conf.gro
    cd $cwd
}

function f02_relax () {
    dir=02.relax
    if test -d $dir; then
	echo existing dir $dir, do nothing.
	exit
    fi
    mkdir $dir
    cd $dir
    ln -s $cwd/template/equi.mol.mdp ./grompp.mdp
    ln -s $cwd/01.pdb2gmx/conf.gro .
    ln -s $cwd/01.pdb2gmx/topol.top .
    gmx grompp &> gmx_grompp.log
    gmx mdrun -v &> gmx_mdrun.log
    cd $cwd
}

function f03_reg_zeta_theta () {
    dir=03.regzt
    if test -d $dir; then
	echo existing dir $dir, do nothing.
	exit
    fi
    mkdir $dir
    cd $dir
    
    ln -s $cwd/template/nvt.mdp ./grompp.mdp
    ln -s $cwd/02.relax/confout.gro ./conf.gro
    ln -s $cwd/02.relax/topol.top ./topol.top 
    
    natoms=`cat conf.gro | head -n 2 | tail -n 1`
    nalpha=`echo "($natoms - 12) / 10" | bc`
    $plm_tool_dir/gen.plumed.py res -s 100 -a $nalpha -i 5 7 9 15 7 9 15 17 6 5 7 9 9 15 17 18 7 9 15 11 7 9 15 10 > plumed.dat

    for ii in `seq 0 $(($nalpha - 1))`;
    do
	pi=`printf %02d $ii`
	angle=$pi-00
	sed -i "/res-$angle/s/AT=.*/AT=-1.5/g" plumed.dat 
	angle=$pi-01
	sed -i "/res-$angle/s/AT=.*/AT=-0.5/g" plumed.dat 
	angle=$pi-04
	sed -i "/res-$angle/s/AT=.*/AT=-2.09/g" plumed.dat 
	angle=$pi-05
	sed -i "/res-$angle/s/AT=.*/AT=2.09/g" plumed.dat 
    done

    gmx grompp &> gmx_grompp.log
    gmx mdrun -v -plumed plumed.dat &> gmx_mdrun.log

    cd $cwd
}

function f04_solvate () {
    dir=04.solvate
    if test -d $dir; then
	echo existing dir $dir, do nothing.
	exit
    fi
    mkdir $dir
    cd $dir
    
    ln -s $cwd/03.regzt/confout.gro ./protein.gro
    gmx solvate -cp protein.gro -cs spc216.gro -o conf.gro -box $box $box $box &> gmx_solvate.log
    nsol=`grep "SOL" gmx_solvate.log  | grep Number | awk '{print $5}'`
    echo number of SOL is $nsol

    ln -s $cwd/template/npt.equi.mdp ./grompp.mdp
    cp -L $cwd/03.regzt/topol.top ./topol.top
    echo "SOL $nsol" >> topol.top

    gmx grompp &> gmx_grompp.log
    gmx mdrun -v &> gmx_mdrun.log

    cd $cwd
}

function f05_equi () {
    dir=05.equi
    if test -d $dir; then
	echo existing dir $dir, do nothing.
	exit
    fi
    mkdir $dir
    cd $dir

    ln -s $cwd/04.solvate/confout.gro ./conf.gro
    ln -s $cwd/04.solvate/topol.top ./topol.top
    ln -s $cwd/template/npt.run.mdp ./grompp.mdp
    
    gmx grompp &> gmx_grompp.log
    gmx mdrun -v &> gmx_mdrun.log

    cd $cwd
}

echo cut pdb 11
cut_pdb_11 $nalpha
echo using box $box
echo doing 01 pdb2gmx
f01_pdb2gmx 
echo doing 02 relax
f02_relax
echo doing 03 reg zeta theta
f03_reg_zeta_theta
echo doing 04 solvate
f04_solvate
echo doing 05 equi
f05_equi
