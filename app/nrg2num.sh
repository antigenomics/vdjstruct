#!/bin/bash

path=generated/pdbcdr3/energy_mats

for file in `find $path -type f -name "total????.xpm"`
do
	mod=${file#$path'/total'}
	name=${mod/'.xpm'/''}
	echo 'Running '$name'...'
	python src/main/python/pdbstat.py $name nrg
done

#./src/main/gromacs/params/script.sh
