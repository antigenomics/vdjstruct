#!/bin/bash

for file in `find generated/pdbcdr3/dist_mats/cdr3+pep -type f -name "*(0).txt"`
do
#	path=${file%'.'}
	mod=${file#'generated/pdbcdr3/dist_mats/cdr3+pep/'}
#	name=${mod/'.pdb'/''}
	name=${mod/'(0).txt'/''}
	echo 'Running '$name'...'
	
	./src/main/gromacs/params/script.sh $name
done
