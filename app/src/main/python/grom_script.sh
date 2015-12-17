#!/bin/bash

function test {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
    	echo
    	"$@" 2>$logpath
        echo "error with $@"
        echo
        exit 0
    fi
    return 0
}

#split='====================================================n'
name=$3
log_path=fails/$name.log
xpm_path=$2
gromacs_path=$1

cd $gromacs_path 

mkdir fails
test gmx pdb2gmx -f $name.pdb -o $name.gro -water spce -missing -ff oplsaa

test gmx editconf -f $name.gro -o $name.gro -c -d 1.0 -bt cubic
#test gmx solvate -cp $name.gro -cs spc216.gro -o $name.gro -p topol.top
#gmx grompp -f params/ions.mdp -c $name.gro -p topol.top -o ions.tpr
#read
#gmx genion -s ions.tpr -o $name.gro -p topol.top -pname NA -nname CL -nn 8
cd -

seqs=$(python enemat.py $*)

cd -

test gmx grompp -f params/minim.mdp -c $name.gro -p topol.top -n index.ndx -o em$name.tpr
test gmx mdrun -v -deffnm em$name -nb cpu
test gmx enemat -groups groups$name.dat -nlevels 10000 -emat $name.xpm -f em$name.edr
echo $seqs >> total*.xpm
mv em$name.edr groups$name.dat trash

cd -

mkdir $xpm_path
mv ${gromacs_path}/*.xpm $xpm_path
echo 'XPM files can be found in '${xpm_path}

rm -f ${gromacs_path}/*

