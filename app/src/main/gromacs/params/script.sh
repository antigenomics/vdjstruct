#!/bin/bash

function test {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
    	echo
    	"$@" 2>$logpath
        echo "error with $@"
        echo
        exit $status
    fi
    return $status
}

#split='====================================================n'
name=$1
logpath=fails/$name.log
bigpath=src/main/gromacs

cd $bigpath 
rm -f * 

test gmx pdb2gmx -f ../../../../truncpdbs/$name.pdb -o $name.gro -water spce -missing -ff oplsaa

test gmx editconf -f $name.gro -o $name.gro -c -d 1.0 -bt cubic
#test gmx solvate -cp $name.gro -cs spc216.gro -o $name.gro -p topol.top
#gmx grompp -f params/ions.mdp -c $name.gro -p topol.top -o ions.tpr
#read
#gmx genion -s ions.tpr -o $name.gro -p topol.top -pname NA -nname CL -nn 8
cd -

seqs=$(python src/main/python/enemat.py $name)
#python src/main/python/enemat.py $name
cd $bigpath
test gmx grompp -f params/minim.mdp -c $name.gro -p topol.top -n index.ndx -o em$name.tpr
test gmx mdrun -v -deffnm em$name -nb cpu
test gmx enemat -groups groups$name.dat -nlevels 10000 -emat $name.xpm -f em$name.edr
echo $seqs >> total*.xpm

mv *.xpm ../../../generated/pdbcdr3/energy_mats
mv em$name.edr groups$name.dat trash
rm -f *

cd -
