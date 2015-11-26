#!/bin/bash

function test {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
    	echo
        echo "error with $@"
        echo
        exit $status
    fi
    return $status
}

for file in `find ../pdbs -type f -name "*.pdb"`
do
	test pdbfixer $file --verbose --add-atoms=heavy --replace-nonstandard --add-residues --output=../fixedpdbs/${file#'../pdbs/'}
done
