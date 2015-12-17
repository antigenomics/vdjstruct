import pdbmod
import pandas as pd
import sys
from Bio.PDB import PDBParser

pdb_path = '../pdbs/'
data_path = '../pdbs/final.annotations.txt'

def Nrg(name):
	print "\nProtein name: " + name + "\n"
	info = data[data['pdb_id'] == name]
	if info.empty:
		return 0
	xpm_path = 'generated/pdbcdr3/energy_mats/' + 'total' + name + '.xpm'
	c = pdbmod.Interaction(info, pdb_path + name + '.pdb')
	if not c.calcNrgMatrices(xpm_path, 'peptide', ('alpha', 'CDR3'), ('beta', 'CDR3')):
		return 0
	if not c.writeInFile_CDR3_Pept_Nrg():
		return 0
	return c

def Dist(name):
	print "\nProtein name: " + name + "\n"
	info = data[data['pdb_id'] == name]
	if info.empty:
		return 0
	c = pdbmod.Interaction(info, pdb_path + name + '.pdb')
	if not c.calcDistMatrices('peptide', ('alpha', 'CDR3')):
		return 0
	if not c.calcDistMatrices('peptide', ('beta', 'CDR3')):
		return 0
	if not c.writeInFile_CDR3_Pept_Dist():
		return 0
	return c
	
def CompressPDB(name):
	print "\nProtein name: " + name + "\n"
	info = data[data['pdb_id'] == name]
	if info.empty:
		return 0
	c = pdbmod.Interaction(info, pdb_path + name + '.pdb')
	if not c.pushToPDB('../truncpdbs/', 'peptide', ('alpha', 'CDR3'), ('beta', 'CDR3')):
		return 0
	return c
	
def Iterate(items, FUN):
	return filter(lambda x: x, [FUN(i) for i in items])
		
	
def Err():
	print 'Wrong flag (nrg/dist/comp)'
	exit(1)

names = map(lambda x: x[-8:-4], sys.argv[1:-1])
flag = sys.argv[-1]

switch = {'nrg': lambda x: Iterate(x, Nrg),
	'dist': lambda x: Iterate(x, Dist),
	'comp': lambda x: Iterate(x, CompressPDB)}
exception = lambda x : Err();
FUN = switch.get(flag, exception)
data = pd.DataFrame(pd.read_table(data_path, sep='\t'))

cpx = FUN(names)
if not cpx:
	print 'ERROR OCCURED!'
#for x in cpx:
#	print x




