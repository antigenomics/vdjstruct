#!/bin/bash

path=../pdbs

mis=0
for file in `find $path -type f -name "*.pdb"`
do
	mod=${file#$path'/'}
	name=${mod/'.pdb'/''}
	echo 'Running '$name'...'
	python src/main/python/pdbstat.py $name dist
	let mis=$mis+$?
done
echo $mis
#./src/main/gromacs/params/script.sh
