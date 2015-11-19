#!/bin/bash

path=generated/pdbcdr3/dist_mats/cdr3+pep

for file in `find $path -type f -name "*(0).txt"`
do
#	path=${file%'.'}
	mod=${file#$path'/'}
#	name=${mod/'.pdb'/''}
	name=${mod/'(0).txt'/''}
	echo 'Running '$name'...'
	
	./src/main/gromacs/params/script.sh $name
done
