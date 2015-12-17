import pdbmod
import pandas as pd
import sys
from Bio.PDB import PDBParser





def Iterate(items, FUN):
	return filter(lambda x: x, [FUN(i) for i in items])
		
	
def Err():
	print 'Wrong flag (nrg/dist/comp)'
	exit(1)

names = map(lambda x: x[-8:-4], sys.argv[1:-1])
flag = sys.argv[-1]

switch = {'nrg': lambda x: Iterate(x, Nrg),
	'dist': lambda x: Iterate(x, Dist),
	'comp': lambda x: Iterate(x, CompressPDB),
	'conv': lambda x: ConvertTable(x)}
exception = lambda x : Err();
FUN = switch.get(flag, exception)
data = pdbmod.parseDataFile("../pdbs/data.txt")

cpx = FUN(names)
