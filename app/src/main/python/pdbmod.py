from Bio.PDB import PDBParser
from Bio.PDB import PPBuilder
from Bio.PDB import Polypeptide

#flag = True if you wanna write useless stuff
flag = False

class Interaction:
	def __init__(self, struct, chains, cdr3_a, cdr3_b):
		self.__name = struct.get_id()		
		self.__struct = struct
		self.__chains = chains
		self.__cdr3_a_seq = cdr3_a
		self.__cdr3_b_seq = cdr3_b
		self.__a_flanked = []
		self.__b_flanked = []
		self.__p_flanked = []
		self.__peptide = []
		self.__cdr3_a = []
		self.__cdr3_b = []
		self.__d_matrix_a = []
		self.__d_matrix_b = []
		self.__e_matrix_a = []
		self.__e_matrix_b = []
		self.__info = ''
		self.getCDR3Indices()
		self.definePeptideChain()
		self.calcDistMatrices()
		self.verbose = False
		
	def __str__(self):
		s = '\nprotein name: ' + self.__name + \
			'\npeptide: ' + self.getPeptideSeq() + \
			'\ncdr3 alpha: ' + self.__cdr3_a_seq + \
			'\ncdr3 beta: ' + self.__cdr3_b_seq + '\n'
		return s
		
	def printerr(self, err):
		if self.verbose:
			print err
		
	def isEmpty(self):
		if not self.__cdr3_a:
			return True
		if not self.__cdr3_b:
			return True
		if not self.__peptide:
			return True
		
	def getCDR3Indices(self):
		ppb=PPBuilder()
		res = []
		bltpep = ppb.build_peptides(self.__struct[0])
		cdr3_a_beg = 1
		cdr3_a_end = 1
		cdr3_b_beg = 1
		cdr3_b_end = 1
		for pp in bltpep: 
			s = str(pp.get_sequence())
			ind = s.find(self.__cdr3_a_seq, 0, len(s))
			if (ind != -1):
				for i in range(ind, ind + len(self.__cdr3_a_seq)):
					res.append(pp[i])
				self.__cdr3_a = res
				cdr3_a_beg = cdr3_a_beg + ind
				cdr3_a_end = cdr3_a_beg + len(self.__cdr3_a_seq) - 1
				self.__a_flanked = [cdr3_a_beg, cdr3_a_end]
				break
			cdr3_a_beg = cdr3_a_beg + len(pp)
				
		if not res:	
			line = '\n' + self.__cdr3_a_seq + ' not found in '+str(self.__struct.get_id()) + '!\n'
			self.__info = self.__info + line
			self.printerr('\ngetCDR3Indices(): ' + line)
		
		res = []
		for pp in bltpep: 
			s = str(pp.get_sequence())
			ind = s.find(self.__cdr3_b_seq, 0, len(s))
			if (ind != -1):
				for i in range(ind, ind+len(self.__cdr3_b_seq)):
					res.append(pp[i])
				self.__cdr3_b = res
				cdr3_b_beg = cdr3_b_beg + ind
				cdr3_b_end = cdr3_b_beg + len(self.__cdr3_b_seq) - 1
				self.__b_flanked = [cdr3_b_beg, cdr3_b_end]
				break
			cdr3_b_beg = cdr3_b_beg + len(pp)
			
		if not res:
			line = '\n'+self.__cdr3_b_seq+' not found in '+str(self.__struct.get_id())+'!\n'
			self.printerr('\ngetCDR3Indices(): ' + line)
			self.__info = self.__info + line
				
	def definePeptideChain(self):
		l = 'INFINITY'
		for chid in self.__chains:
			buf = len(PPBuilder().build_peptides(self.__struct[0][chid])[0])
			if (buf <= l):
				l = buf
				i = chid
		pp = list(self.__struct[0][i])
		if (len(pp) > 20):
			line = self.__name + '\t: TOO MANY AMINO ACIDS (' + str(len(pp)) + ') TO BE A PEPTIDE :(\n'
			self.printerr('\ndefinePeptideChain(): ' + line)
			self.__info = self.__info + line
		else:
			rnum = 1
			for chain in self.__struct[0]:
				if chain != self.__struct[0][i]:
					rnum = rnum + len(PPBuilder().build_peptides(chain)[0])
				else:
					break
			cdr3_p_beg = rnum	
			newpp = []
			for r in pp:
				if (Polypeptide.is_aa(r.get_resname(), standard=True)):
					newpp.append(r)
			self.__peptide = newpp
			cdr3_p_end = cdr3_p_beg + len(self.__peptide) - 1
			self.__p_flanked = [cdr3_p_beg, cdr3_p_end]
						
	def calcDistMatrices(self):
		mat = []
		if not self.__peptide:
			self.printerr('\ncalcDistMatrices(): PEPTIDE IS EMPTY\n')
			return
		if not self.__cdr3_a:
			self.printerr('\ncalcDistMatrices(): ALPHA CDR3 IS EMPTY\n')
		else:
			for res1 in list(self.__peptide):
				mat.append([])	
				for res2 in self.__cdr3_a:	
					mat[len(mat)-1].append(residuesMinDist(res1, res2))
			self.__d_matrix_a = mat
		mat = []
		if not self.__cdr3_b:
			self.printerr('\ncalcDistMatrices(): BETA CDR3 IS EMPTY\n')
		else:
			for res1 in list(self.__peptide):
				mat.append([])		
				for res2 in self.__cdr3_b:
					mat[len(mat)-1].append(residuesMinDist(res1, res2))
			self.__d_matrix_b = mat
			
	def getNrgMatrices(self, path):
		# check
		a, b, annot = extractXPM(path)
		slist = [self.__name, self.getPeptideSeq(), self.getCDR3AlphaSeq(), self.getCDR3BetaSeq()]
		if annot != slist:
			line = slist + '\n' + annot +  '\nMISMATCH\n'
			self.printerr('\nextractXPM(): \n' + line)
			self.__info = self.__info + line
			return

		self.__e_matrix_a = a
		self.__e_matrix_b = b
		
	def getPeptideSeq(self):
		if not self.__peptide:
			self.printerr('\ngetPeptideSeq(): PEPTIDE (' + self.__name +') IS EMPTY\n')
			return	
		s = ''
		for r in list(self.__peptide):
			s = s + Polypeptide.three_to_one(r.get_resname())
		return s
		
	def getCDR3AlphaSeq(self):
		return self.__cdr3_a_seq
		
	def getCDR3BetaSeq(self):
		return self.__cdr3_b_seq
	
	@staticmethod	
	def matToStr(pep, cdr3, mat):
		s = '/'
		for r in cdr3:
			s = s + '\t' + Polypeptide.three_to_one(r.get_resname())
		s = s + '\n'
		for i in range(0, len(mat)):
			s = s + Polypeptide.three_to_one(pep[i].get_resname())
			for j in range(0, len(mat[i])):
				s = s +'\t' + str(mat[i][j])
			s = s + '\n'
		return s
		
	def a_matToStr(self):
		if not self.__d_matrix_a:
			self.printerr('\na_matToStr(): ALPHA MATRIX IS EMPTY\n')
			return
		s = Interaction.matToStr(self.__peptide, self.__cdr3_a, self.__d_matrix_a)
		return s
		
	def b_matToStr(self):
		if not self.__d_matrix_b:
			self.printerr('\nb_matToStr(): BETA MATRIX IS EMPTY\n')
			return
		s = Interaction.matToStr(self.__peptide, self.__cdr3_b, self.__d_matrix_b)
		return s
		
	def a_e_matToStr(self):
		if not self.__e_matrix_a:
			self.printerr('\na_matToStr(): ALPHA MATRIX IS EMPTY\n')
			return
		s = Interaction.matToStr(self.__peptide, self.__cdr3_a, self.__e_matrix_a)
		return s
		
	def b_e_matToStr(self):
		if not self.__e_matrix_b:
			self.printerr('\nb_matToStr(): BETA MATRIX IS EMPTY\n')
			return
		s = Interaction.matToStr(self.__peptide, self.__cdr3_b, self.__e_matrix_b)
		return s
	
	def writeInFile_CDR3_Pept(self, f):
		if self.isEmpty():
			f.write(self.__info)
			f.write('\n')
			return
			
		for i in range(0, 2):
			f2 = open('generated/pdbcdr3/dist_mats/cdr3+pep/'+self.__name+'('+str(i)+').txt', 'w')
			if (i == 0):
				f.write(self.a_matToStr())
				f2.write(self.a_matToStr())
			if (i == 1):
				f.write(self.b_matToStr())
				f2.write(self.b_matToStr())
			f.write('\n')
			f2.write('\n')
			f2.close()
			
	def writeInFile_CDR3_Pept_Nrg(self):
		if self.isEmpty():
			self.printerr(self.__info)
			return
			
		for i in range(0, 2):
			f = open('generated/pdbcdr3/energy_mats/'+self.__name+'('+str(i)+').txt', 'w')
			if (i == 0):
				f.write(self.a_e_matToStr())
			if (i == 1):
				f.write(self.b_e_matToStr())
			f.write('\n')
			f.close()
			
	def indicesToStr(self):
		s = ''
		if self.isEmpty():
			f.write(self.__info)
			return s
		for r in self.__peptide:
			s = s + str(r.get_id())+' '
		s = s + '\n'
		for r in self.__cdr3_a:
			s = s + str(r.get_id())+' '
		s = s + '\n'
		for r in self.__cdr3_b:
			s = s + str(r.get_id())+' '
		s = s + '\n\n'
		return s
		
	def printIndices(self, f):
		f.write(self.indicesToStr())
		
	def utilGromacs(self):
		print [self.__p_flanked, self.__a_flanked, self.__b_flanked]
		
	
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
		
def residuesCaDist(res2, res1):
	if (res1 == res2):
		return 0;
	else:
		return res1['CA']-res2['CA']
		

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
	
def extractXPM(path):
	xpm = open(path, 'r')
	lines = list(xpm)
	xpm.close()
	ind = 0
	for line in lines:
		pline = line.split()
		if pline[0] == 'static':
		    break
		ind += 1

	ind += 1 
	print lines[ind].split('"')[0]
	snum = int(lines[ind].split('"')[1].split()[3])
	ind += 1 

	alphabet = {}
	while(True):
		if lines[ind][0] == '/':
		    break
		pline = lines[ind].split('"')
		code = pline[1].split()[0]
		alphabet.update([(code, pline[3])])
		ind += 1
	ind += 2

	mat = []
	while(True):
		if lines[ind][0] == '/':
		    break
		line = lines[ind].split('"')[1]
		l = []
		ins = ''
		for i in range(0, len(line)):
		    ins += line[i]
		    if (i+1)%snum == 0:
			l.append(float(alphabet[ins]))
			ins = ''
		mat.append(l)
		ind += 1
		
	lline = lines[ind]
	annot = lline.replace('/', ' ').replace('*', ' ').split()
	mat = list(reversed(mat))
		
	a_mat = []
	a_bs = [len(annot[1]), len(annot[1]) + len(annot[2])]
	for line in mat[:len(annot[1])]:
		a_mat.append(line[a_bs[0]:a_bs[1]])
		
	b_mat = []
	b_bs = [a_bs[1], a_bs[1] + len(annot[3])]
	for line in mat[:len(annot[1])]:
		b_mat.append(line[b_bs[0]:b_bs[1]])
		
	return a_mat, b_mat, annot

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
