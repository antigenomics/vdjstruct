from Bio.PDB import PDBParser
from Bio.PDB import PPBuilder
from Bio.PDB import Polypeptide

#flag = True if you wanna write useless stuff
flag = False

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
		
def residuesCaDist(res1, res2):
	dist = 'INFINITY'
	if (res1 == res2):
		return 0;
	else:
		return res1['CA']-res2['CA']

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

def calcDistMatrix(pname, pept1, pept2, symflag):
	mat = [[pname]]
	indlist = []
	if (not symflag):
		#mat.append(['/'])
		mat[0] = ['/']
	else:
		#mat.append([])
		mat[0] = []
		
	if (len(pept1) == 0):
		if (flag):
			mat[len(mat)-1].append('CDR3 NOT FOUND!')
		return mat
	for res1 in pept1:
		name = res1.get_resname()
		if (Polypeptide.is_aa(name, standard=True)):
			mat[len(mat)-1].append(Polypeptide.three_to_one(name))
	for res2 in pept2:		
		name = res2.get_resname()
		if (Polypeptide.is_aa(name, standard=True)):
			mat.append([])
			if (not symflag):
				mat[len(mat)-1].append(Polypeptide.three_to_one(name))
			for res1 in pept1:
				if (Polypeptide.is_aa(res1.get_resname())):
					mat[len(mat)-1].append(residuesCaDist(res1, res2))
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
		for j in range(0, len(m[i])):
			f.write(str(m[i][j])+'\t')
		f.write('\n')

def writeInFile_CDR3_Pept(parser, datadist, item, f):
	structure = parser.get_structure(item, '../pdbs/'+item+'.pdb')
	cdr3seqlist = datadist[item]
	peptide = structure[0][definePeptideChain(cdr3seqlist[0], structure)]
	if (len(peptide) > 20):
		if (flag == True):
			line = item+'\tTOO MANY AMINO ACIDS ('+str(len(peptide))+') TO BE A PEPTIDE :(\n\n'
			print line
			f.write(line)
	else:
		cdr3reslist = []
		for i in range(1, len(cdr3seqlist)):
			cdr3reslist.append(getCDR3Indices(structure, cdr3seqlist[i]))
		i = 0
		for cdr3res in cdr3reslist:
			f2 = open('generated/pdbcdr3/dist_mats/cdr3+pep/'+item+'('+str(i)+').txt', 'w')
			m = calcDistMatrix(item, cdr3res, list(peptide), False)
			fwriteMatrix(f2, m)
			f2.write('\n')
			f2.close()
			f.write(item+'\n')	
			fwriteMatrix(f, m)
			f.write('\n')
			i += 1
			
def writeInFile_CDR3_CDR3(parser, datadist, item, f):
	structure = parser.get_structure(item, '../pdbs/'+item+'.pdb')
	cdr3seqlist = datadist[item]
	cdr3reslist = []
	for i in range(1, len(cdr3seqlist)):
		cdr3reslist.append(getCDR3Indices(structure, cdr3seqlist[i]))
	i = 0
	for cdr3res in cdr3reslist:
		f2 = open('generated/pdbcdr3/dist_mats/cdr3+cdr3/'+item+'('+str(i)+').txt', 'w')
		m = calcDistMatrix(item, cdr3res, cdr3res, True)
		fwriteMatrix(f2, m)
		f2.write('\n')
		f2.close()
		f.write(item+'\n')
		fwriteMatrix(f, m)
		f.write('\n')
		i += 1


