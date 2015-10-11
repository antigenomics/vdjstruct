import pdbmod
from Bio.PDB import PDBParser

parser = PDBParser()
datadist = pdbmod.parseDataFile("../pdbs/data.txt")
rts = open('generated/pdbcdr3/cdr3_pep_dist_mats.txt', 'w')
for item in datadist:	
	pdbmod.writeInFile_CDR3_Pept(parser, datadist, item, rts)
rts.close()
