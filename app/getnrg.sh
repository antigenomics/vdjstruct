#!/bin/bash

path=generated/pdbcdr3/dist_mats/cdr3+pep
path=src/main/gromacs/fails

function test {
    "$@"
    local status=$?
    if [ $status == 0 ]; then
    	rm src/main/gromacs/fails/$name.log
        continue
    fi
    return $status
}

for file in `find $path -type f -name "*3vxr.log"`  #`find $path -type f -name "*(0).txt"`
do
	#mod=${file#$path'/'}
	#name=${mod/'(0).txt'/''}
	mod=${file#$path'/'}
	name=${mod/'.log'/''}
	echo 'Running '$name'...'
	
	test ./src/main/gromacs/params/script.sh $name
	
	echo
	echo "Error occured: processing $name..."
	echo
	
	read
	cd ../fixedpdbs
	python fix.py $name
	err=$?
	if [ "$err" == '1' ]; then
		echo "Seems there are problems with .edr file $name..."
		echo
		cp ../pdbs/$name.pdb ./
	fi 
	cd -
	read
	test ./src/main/gromacs/params/script.sh $name
	
	echo
	echo "Error occured, while getting energies"
done
