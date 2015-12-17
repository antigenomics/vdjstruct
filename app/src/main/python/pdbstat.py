import pdbmod
import pandas as pd
import sys
from Bio.PDB import PDBParser

#==================================================================================
# For calculating energies grom_script.sh must be in the same folder with this file
#
#			    
# Usage example: pdbstat.py name1.pdb name2.pdb name3.pdb table
#
# this will create summary table in the folder specified by result_table_dir
# with the regions listed in the variables main_region and regions
#==================================================================================
# variants for final flag: nrg, dist, comp, table
#==================================================================================

pdb_path = '../../../../pdbs'				# where .pdb files lie
data_path = '../../../../pdbs/final.annotations.txt'	# annotation file
nrg_dir = '../../../generated/pdbcdr3/energy_mats' 	# where to write nrg matrices
dist_dir = '../../../generated/pdbcdr3/dist_mats'	# where to write dist matrices
gromacs_dir = '../gromacs'				# sandbox for gromacs actions
xpm_dir = '../gromacs/xpm'				# folder with .xmp files (result of gromacs activity)
result_table_dir = '.'					# where to write summary table

# specify regions you want to calculate intercation between (for Table(name))
main_region = 'peptide' # other variants of region names : ('alpha', 'CDR3') or ('alpha', 'CDR1')
			# or ('beta', 'CDR1') etc. according to the names in the final.annotations table
			
	    # first region      second region     enter as much region as you wish and launch
regions = [('alpha', 'CDR3'), ('beta', 'CDR3')] 


def Nrg(name): # store a couple of energy matrices for peptide|CDR3alpha (0) and peptide|CDR3beta (1) to nrg_dir
	print "\nProtein name: " + name + "\n"
	info = data[data['pdb_id'] == name]
	if info.empty:
		print name, ' WAS NOT FOUND IN THE TABLE'
		return 0
	c = pdbmod.Interaction(info, pdb_path + '/' + name + '.pdb')
	if not c.calcNrgMatrices(gromacs_dir, xpm_dir, 'peptide', ('alpha', 'CDR3'), ('beta', 'CDR3')):
		return 0
	if not c.writeInFile_CDR3_Pept_Nrg(nrg_dir):
		return 0
	return c

def Dist(name):	# store a couple of distance matrices for peptide|CDR3alpha (0) and peptide|CDR3beta (1) to dist_dir
	print "\nProtein name: " + name + "\n"
	info = data[data['pdb_id'] == name]
	if info.empty:
		print name, 'WAS NOT FOUND IN THE TABLE'
		return 0
	try:
		c = pdbmod.Interaction(info, pdb_path + '/' + name + '.pdb')
	except IOError:
		print 'PDB FILE ', pdb_path, 'DOES NOT EXIT'
		return 0
	if not c.calcDistMatrices('peptide', ('alpha', 'CDR3')):
		return 0
	if not c.calcDistMatrices('peptide', ('beta', 'CDR3')):
		return 0
	if not c.writeInFile_CDR3_Pept_Dist(dist_dir):
		return 0
	return c
	
def CompressPDB(name): # store truncated pdbs contaning only specified regions to ../truncpdbs
	print "\nProtein name: " + name + "\n"
	info = data[data['pdb_id'] == name]
	if info.empty:
		print name, 'WAS NOT FOUND IN THE TABLE'
		return 0
	c = pdbmod.Interaction(info, pdb_path + '/' + name + '.pdb')
	if not c.pushToPDB('../../../../truncpdbs', 'peptide', ('alpha', 'CDR3'), ('beta', 'CDR3')):
		return 0
	return c
	
def Table(name): # create big table with pairwise energy and distance
	print "\nProtein name: " + name + "\n"
	info = data[data['pdb_id'] == name]
	if info.empty:
		print name, 'WAS NOT FOUND IN THE TABLE'
		return 0
	try:
		c = pdbmod.Interaction(info, pdb_path + '/' + name + '.pdb')
	except IOError:
		print 'PDB FILE ', pdb_path, 'DOES NOT EXIT'
		return 0
		
	for region in regions:
		if not c.calcDistMatrices(main_region, region):
			return 0
	if not c.calcNrgMatrices(gromacs_dir, xpm_dir, main_region, *regions):
		return 0
	if not c.overallTable(main_region, *regions):
		return 0
	return c
	
def Iterate(items, FUN):	# run one of upper functions for several proteins
	return filter(lambda x: x != 0, [FUN(i) for i in items])
		
	
def Err():
	print 'Wrong flag (nrg/dist/comp/table)'
	exit(1)


names = map(lambda x: x[-8:-4], sys.argv[1:-1])
flag = sys.argv[-1]

switch = {'nrg': lambda x: Iterate(x, Nrg),
	'dist': lambda x: Iterate(x, Dist),
	'comp': lambda x: Iterate(x, CompressPDB), 
	'table': lambda x: Iterate(x, Table)}
exception = lambda x : Err();
FUN = switch.get(flag, exception)
data = pd.DataFrame(pd.read_table(data_path, sep='\t'))
if len(sys.argv) == 3:
	if sys.argv[1] == 'all':
		names = list(set(data['pdb_id'].tolist()))

cpx = FUN(names)
if not cpx:
	print 'ERROR OCCURED!'
	
try:
	if flag == 'table':
		result_table = pd.concat([x.summary_table for x in cpx])
		result_table.to_html(result_table_dir + '/summary_table.html')
except ValueError:
	print 'NOTHING TO RETURN'




