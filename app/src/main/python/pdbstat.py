import pdbmod
import sys
from Bio.PDB import PDBParser

def Nrg(name):
	info = pdbmod.getProtein(data, name)
	structure = PDBParser().get_structure(name, '../pdbs/'+name+'.pdb')
	c = pdbmod.Interaction(structure, *info)
	c.getNrgMatrices('generated/pdbcdr3/energy_mats/' + 'total' + name + '.xpm')
	c.writeInFile_CDR3_Pept_Nrg()
	return c

def Dist(name):
	info = pdbmod.getProtein(data, name)
	structure = PDBParser().get_structure(name, '../pdbs/'+name+'.pdb')
	c = pdbmod.Interaction(structure, *info)
	c.calcDistMatrices()
	c.writeInFile_CDR3_Pept_Dist()
	return c
	
def CompressPDB(name):
	info = pdbmod.getProtein(data, name)
	structure = PDBParser().get_structure(name, '../pdbs/'+name+'.pdb')
	c = pdbmod.Interaction(structure, *info)
	c.pushToPDB()
	return c
	
def Err():
	print 'Wrong flag (nrg/dist/comp)'
	exit(1)

item = sys.argv[1]
flag = sys.argv[2]

switch = {'nrg': lambda x: Nrg(x),
	'dist': lambda x: Dist(x),
	'comp': lambda x: CompressPDB(x)}
exception = lambda x : Err();

print "Protein name: " + item + "\n"
FUN = switch.get(flag, exception)

data = pdbmod.parseDataFile("../pdbs/data.txt")
	
cpx = FUN(item)

#cpx.convertTable()

