#!/bin/bash

function test {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
        echo "error with $@" > $logpath
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

test gmx pdb2gmx -f ../../../../pdbs/$name.pdb -o $name.gro -water spce -ff oplsaa -missing
test gmx editconf -f $name.gro -o $name.gro -c -d 1.0 -bt cubic
test gmx solvate -cp $name.gro -cs spc216.gro -o $name.gro -p topol.top

#gmx grompp -f params/ions.mdp -c $name.gro -p topol.top -o ions.tpr
#read
#gmx genion -s ions.tpr -o $name.gro -p topol.top -pname NA -nname CL -nn 8
#read
cd -

python src/main/python/enemat.py $name

cd $bigpath
test gmx grompp -f params/minim.mdp -c $name.gro -p topol.top -n index.ndx -o em$name.tpr
test gmx mdrun -v -deffnm em$name
test gmx enemat -groups groups$name.dat -nlevels 200 -emat $name.xpm -f em$name.edr

mv *.xpm ../../../generated/pdbcdr3/energy_mats
mv em$name.edr groups$name.dat trash
rm -f *

cd -
