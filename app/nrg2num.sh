#!/bin/bash

path=generated/pdbcdr3/energy_mats

python src/main/python/pdbstat.py `find $path -type f -name "total????.xpm"` nrg

#./src/main/gromacs/params/script.sh
