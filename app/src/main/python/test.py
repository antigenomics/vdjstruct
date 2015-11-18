from Bio.PDB import PDBParser
from Bio.PDB import PPBuilder
from Bio.PDB import Polypeptide

item = '2bnr'	
structure = PDBParser().get_structure(item, '../pdbs/'+item+'.pdb')
ppb=PPBuilder()
peps = ppb.build_peptides(structure)

print structure.get_id()
print peps[0]
#print peps[0][1:-3]
print peps[0][3:9]
p = peps[0][3:9]
print peps[0][1].get_resname()
