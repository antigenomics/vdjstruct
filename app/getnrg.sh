#!/bin/bash

#path=generated/pdbcdr3/dist_mats/cdr3+pep
#path=src/main/gromacs/fails
path=../pdbs
logpath=getnrgfails/
function test {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
    	echo
    	"$@" >> getnrgfails/$name.log
        echo "error with $@"
        echo
        continue
        #exit $status
    fi
    return $status
}

for file in `find $path -type f -name "*.pdb"`  #`find $path -type f -name "*(0).txt"`
do
	mod=${file#$path'/'}
	name=${mod/'.pdb'/''}
	echo 'Running '$name'...'
	
	test python src/main/python/pdbstat.py $name comp
        ./src/main/gromacs/params/script.sh $name
done
