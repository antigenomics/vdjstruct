
# coding: utf-8

# In[237]:

from Bio.PDB import PDBParser
from Bio.PDB import Polypeptide
from gromacs import cbook
import gromacs.cbook
import pdbmod
import sys


# In[238]:

def parseGro(f):
    
    fp = open(f, 'r')
    pp = list(fp)
    res = ''
    atlist = []
    aas = []
    prevnum = '0'
    for i in range(2, len(pp)):
        line = pp[i].split()
        num = line[0][:-3]
        aa = line[0][-3:]
        
        if not Polypeptide.is_aa(aa):
            break;
            
        if (num != prevnum):
            shortaa = Polypeptide.three_to_one(aa)
            atlist.append(i - 1)
            aas.append(aa + '_' + num)
            prevnum = num
            res = res + shortaa
            
    fp.close()
    return res, atlist, aas


# In[239]:

def findSeqsInGro(seq, pep, cdr3_a, cdr3_b):
    pbeg = seq.find(pep) + 1
    pend = pbeg + len(pep) - 1
    abeg = seq.find(cdr3_a) + 1
    aend = abeg + len(cdr3_a) - 1
    bbeg = seq.find(cdr3_b) + 1
    bend = bbeg + len(cdr3_b) - 1
    return [pbeg, pend], [abeg, aend], [bbeg, bend]


# In[240]:

def appendSeqs(path, idlist, atlist, aas):
    fl = open(path, 'a')
    groups = gromacs.cbook.get_ndx_groups(path)
    names = []
    for pair in idlist:
        for aaid in range(pair[0] - 1, pair[1]):
            name = 'r_' + str(pair[0]) + '-' + str(pair[1]) + '_' + aas[aaid]
            names.append(name)
            if any(group['name'] == name for group in groups):
                continue
            fl.write('[ ' + name + ' ]\n')
            for item in range(atlist[aaid], atlist[aaid + 1]):
                fl.write(str(item) + ' ')
            fl.write('\n')
    fl.close()
    return names


# In[241]:

def appendToMDP(path, grnames):
    replaceText = 'energygrps\t='
    for name in grnames:
        replaceText = replaceText + ' ' + name
    replaceText = replaceText + '\n'
    
    fl = open(path, 'r')
    text = fl.read()
    fl.close()
    
    fl = open(path, 'r')
    s = ''
    for line in fl:
        if line.split()[0] == 'energygrps':
            s = line
            break
    fl.close()
    
    fl = open(path, 'w')
    if not s:
        fl.write(text + replaceText)
    else:
        fl.write(text.replace(s, replaceText))
    fl.close()


# In[242]:

def appendToGroupsDat(path, names):
    fl = open(path, 'w')
    fl.write(str(len(names)) + '\n')
    for name in names:
        fl.write(name + '\n')


# In[243]:

def appendSeqs(path, idlist, atlist, aas):
    fl = open(path, 'a')
    groups = gromacs.cbook.get_ndx_groups(path)
    names = []
    for pair in idlist:
        for aa in range(pair[0] - 1, pair[1]):
            name = 'r_' + str(pair[0]) + '-' + str(pair[1]) + '_' + aas[aa]
            names.append(name)
            if any(group['name'] == name for group in groups):
                continue
            fl.write('[ ' + name + ' ]\n')
            for atom in range(atlist[aa], atlist[aa + 1]):
                fl.write(str(atom) + ' ')
            fl.write('\n')
    fl.close()
    return names

item = sys.argv[1]

indpath = 'src/main/gromacs/'
path = indpath + 'params/'
datadist = pdbmod.parseDataFile("../pdbs/data.txt")
G = gromacs.cbook.IndexBuilder(indpath + item + '.gro')
G.cat(indpath + 'index.ndx')

l = datadist[item]
structure = PDBParser().get_structure(item, '../fixedpdbs/'+item+'.pdb')
protein = pdbmod.Interaction(structure, l[0], l[1], l[2])

st, atlist, aas = parseGro(indpath + item + '.gro')

pseq = protein.getPeptideSeq()
print pseq
aseq = protein.getCDR3AlphaSeq()
print aseq
bseq = protein.getCDR3BetaSeq()
print bseq

p, a, b = findSeqsInGro(st, pseq, aseq, bseq)
print p, a, b
names = appendSeqs(indpath + 'index.ndx', [p, a, b], atlist, aas)
appendToMDP(path + 'minim.mdp', names)
appendToGroupsDat(indpath + 'groups' + item + '.dat', names)

print "/*%s*/\n/*%s*/\n/*%s*/\n/*%s*/\n" % (item, pseq, aseq, bseq)


