import pdbmod
from Bio.PDB import PDBParser

datadist = pdbmod2.parseDataFile("../pdbs/data.txt")
rts = open('generated/pdbcdr3/cdr3_pep_dist_mats.txt', 'w')
#
#item = '2bnr'

for item in datadist:	
	l = datadist[item]
	structure = PDBParser().get_structure(item, '../pdbs/'+item+'.pdb')
	c = pdbmod2.Interaction(structure, l[0], l[1], l[2])
c.writeInFile_CDR3_Pept(rts)
#	pdbmod.writeInFile_CDR3_Pept(parser, datadist, item, rts)
rts.close()
print c

c.utilGromacs()
