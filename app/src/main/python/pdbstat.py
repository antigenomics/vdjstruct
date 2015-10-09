from Bio.PDB import PDBParser
from Bio.PDB import PPBuilder
from Bio.PDB import Polypeptide

def residuesMeanDist(res1, res2):
	dist = 0
	num = 0
	if (res1 == res2):
		return 0;
	else:
		for at1 in res1:
			for at2 in res2:
				dist += (at1-at2)
				num += 1
		return dist/num
		
def residuesMinDist(res1, res2):
	dist = 'INFINITY'
	if (res1 == res2):
		return 0;
	else:
		for at1 in res1:
			for at2 in res2:
				buf = at1-at2
				if (buf < dist):
					dist = buf		
		return dist

def getCDR3Indices(struct, cdr3):
	ppb=PPBuilder()
	chid = 0
	res = []
	for pp in ppb.build_peptides(struct): 
		s = str(pp.get_sequence())
		ind = s.find(cdr3, 0, len(s))
		if (ind != -1):
			for i in range(ind, ind+len(cdr3)):
				res.append(pp[i])
			return res
		chid += 1
	print '\n'+cdr3+' not found in '+str(struct.get_id())+'!\n'
	return []

def calcDistMatrix(pname, cdr3, peptide_chain):
	mat = [[]]
	indlist = []
	mat[0].append(pname)
	if (len(cdr3) == 0):
		mat[0].append('CDR3 NOT FOUND!')	
	else:
		for i in range(0, len(cdr3)):
			name = cdr3[i].get_resname()
			if (Polypeptide.is_aa(name)):
				mat[0].append(name)
				indlist.append(i)
		for presidue in peptide_chain:
			name = presidue.get_resname()
			if (Polypeptide.is_aa(name)):
				mat.append([])
				mat[len(mat)-1].append(name)
				for i in indlist:
					mat[len(mat)-1].append(residuesMinDist(cdr3[i], presidue))
	return mat

def definePeptideChain(chains, struct):
	l = len(struct[0][chains[0]])
	for chid in chains:
		buf = len(struct[0][chid])
		if (buf <= l):
			l = buf
			i = chid
	return i

def deleteNonAmino(p):
	for r in p:
		if (Polypeptide.is_aa(r.get_resname())):
			print ""
		else:
			del r

def getPeptideAminoList(p):
	l = []
	for r in p:
		l.append(r.get_resname())
	return l

def parseDataFile(path):
	f = open(path, 'r')
	l = list(f)
	d = {}
	for i in range (1, len(l)):
		s = l[i].split("\t")
		if (len(s) == 20):
			pname = s[0].lower()
			if (d.get(pname) == None):
				d.update({pname : [[s[1]]]})
			else:
				d[pname][0].append(s[1])
			if (s[17].isalpha()):
				d[pname].append(s[17])
	f.close()
	return d
	
def fwriteMatrix(f, m):
	for i in range(0, len(m)):
		for j in range(0, len(m[0])):
			f.write(str(m[i][j])+'\t')
		f.write('\n')

# Main

parser = PDBParser()
datadist = parseDataFile("../pdbs/data.txt")
rts = open('generated/pdbcdr3_proc/dist_mats.txt', 'w')

for item in datadist:	
	structure = parser.get_structure(item, '../pdbs/'+item+'.pdb')
	cdr3seqlist = datadist[item]
	peptide = structure[0][definePeptideChain(cdr3seqlist[0], structure)]
	if (len(peptide) > 20):
		line = item+'\tTOO MANY AMINO ACIDS ('+str(len(peptide))+') TO BE A PEPTIDE :(\n\n'
		print line
		rts.write(line)
	else:
		cdr3reslist = []
		for i in range(1, len(cdr3seqlist)):
			cdr3reslist.append(getCDR3Indices(structure, cdr3seqlist[i]))
		matrices = []
		for cdr3res in cdr3reslist:
			matrices.append(calcDistMatrix(item, cdr3res, peptide))
		for matrix in matrices:
			fwriteMatrix(rts, matrix)
			rts.write('\n')

#builder = PPBuilder()	
#peps = builder.build_peptides(structure[0]['E'])
#for p in peps:
#	print p.get_sequence()
#	print "=========="

rts.close()

