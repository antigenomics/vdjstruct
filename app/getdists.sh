#!/bin/bash

path=../pdbs

python src/main/python/pdbstat.py `find $path -type f -name "*.pdb"` dist

#./src/main/gromacs/params/script.sh
