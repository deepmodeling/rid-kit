#!/bin/bash

GMX=gmx
threads=4

PREP="$GMX grompp"
EXEC="$GMX mdrun -plumed plumed.res.dat -notunepme -nt $threads"
EXEC_CONT="$GMX mdrun -plumed plumed.res.dat -notunepme -nt $threads -cpi state.cpt"

# do equilibrium 
if test ! -f equi.gro; then    
    echo do equi
    $PREP -f equi.mdp -o equi.tpr -po equiout.mdp &> equi.log 
    $EXEC -s equi.tpr -c equi.gro &>> equi.log 
    if test $? -ne 0; then
	echo equi failed, exit
	exit
    fi
    rm -fr state*.cpt md.log ener.edr plm.res.out
fi

if test -f state.cpt; then
    echo prod run from prev
    $EXEC_CONT &> mdrun.log
else
    echo prod run
    $PREP -c equi.gro &> grompp.log
    $EXEC &> mdrun.log
fi    

