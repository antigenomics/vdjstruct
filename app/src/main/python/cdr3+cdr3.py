import pdbmod
from Bio.PDB import PDBParser

parser = PDBParser()
datadist = pdbmod.parseDataFile("../pdbs/data.txt")
rts = open('generated/pdbcdr3/cdr3_cdr3_dist_mats.txt', 'w')
for item in datadist:	
	pdbmod.writeInFile_CDR3_CDR3(parser, datadist, item, rts)
rts.close()
