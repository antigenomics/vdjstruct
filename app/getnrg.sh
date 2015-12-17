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
    	"$@" >> getnrgfails/log
        echo "error with $@"
        echo
        #exit $status
    fi
    return $status
}

test python src/main/python/pdbstat.py `find $path -type f -name "*.pdb"` comp
echo finish
read
for file in `find $path -type f -name "*.pdb"`
do
	mod=${file#$path'/'}
	name=${mod/'.pdb'/''}
	echo 'Running '$name'...'
        ./src/main/gromacs/params/script.sh $name
done
