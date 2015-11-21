#!/bin/bash

function test {
    "$@" 2>$logpath
    local status=$?
    if [ $status -ne 0 ]; then
    	echo
        echo "error with $@"
        echo
        exit $status
    fi
    return $status
}

function test_pdb2gmx {
    "$@" 2>$logpath
    local status=$?
    if [ $status -ne 0 ]; then
    	echo
        echo "error with $@ ------> "
        echo "------> creating primitive topology"
        echo
        test gmx x2top -f ../../../../pdbs/$name.pdb -o topol.top -ff oplsaa 
        flag=0
        #exit $status
    fi
    return $status
}

#split='====================================================n'
name=$1
logpath=fails/$name.log
bigpath=src/main/gromacs
flag=1

cd $bigpath
rm -f * 

test_pdb2gmx gmx pdb2gmx -f ../../../../pdbs/$name.pdb -o $name.gro -water spce -missing -ff oplsaa

if [ $flag == 1 ]; then
	test gmx editconf -f $name.gro -o $name.gro -c -d 1.0 -bt cubic
	test gmx solvate -cp $name.gro -cs spc216.gro -o $name.gro -p topol.top
fi

#gmx grompp -f params/ions.mdp -c $name.gro -p topol.top -o ions.tpr
#read
#gmx genion -s ions.tpr -o $name.gro -p topol.top -pname NA -nname CL -nn 8
#read
cd -

seqs=$(python src/main/python/enemat.py $name)

cd $bigpath
test gmx grompp -f params/minim.mdp -c $name.gro -p topol.top -n index.ndx -o em$name.tpr
test gmx mdrun -v -deffnm em$name
test gmx enemat -groups groups$name.dat -nlevels 10000 -emat $name.xpm -f em$name.edr
echo $seqs >> total*.xpm

mv *.xpm ../../../generated/pdbcdr3/energy_mats
mv em$name.edr groups$name.dat trash
rm -f *

cd -
