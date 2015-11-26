import pdbmod
import sys
from Bio.PDB import PDBParser

datadist = pdbmod.parseDataFile("../pdbs/data.txt")
#rts = open('generated/pdbcdr3/cdr3_pep_dist_mats.txt', 'w')

item = sys.argv[1]
print item 

#for item in datadist:	
l = datadist[item]
structure = PDBParser().get_structure(item, '../pdbs/'+item+'.pdb')
c = pdbmod.Interaction(structure, l[0], l[1], l[2])
#c.writeInFile_CDR3_Pept(rts)
c.getNrgMatrices('generated/pdbcdr3/energy_mats/' + 'total' + item + '.xpm')
#rts.close()
c.writeInFile_CDR3_Pept_Nrg()

#print c.a_e_matToStr()

