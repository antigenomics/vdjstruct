from Bio import SeqIO
from Bio import SeqFeature

def processRec(rec):
	l = rec.split(";")
	score = 0;
	for i in range(0, len(l)):
		l[i] = l[i].split("|")
		newscore = float(l[i][len(l[i])-1])
		if (score <= newscore):
			score = newscore
			result = list([l[i][3], l[i][4], i])
	return result
		
def processLine(string):
	l = string.split("\t")
	result = [l[0]]
	for i in range(6, 10):
		if (l[i] == ""):
			result.append("-1")
			result.append("-1")
			result.append("-1")
		else:
			buf = processRec(l[i])
			result.append(buf[0])
			result.append(buf[1])
			buf1 = l[i-4].split(",")
			buf2 = buf1[buf[2]].split("(")
			result.append(buf2[0])
	return result
	
def writeIn(f, l):
	for i in l:
		f.write(i+"/")
		

align = open("generated/beta_alignments.txt", "r")
cdr = open("generated/beta_cdr3.txt", "r")
fasta = open("generated/beta_filtered.fasta", "r")
result = open("generated/beta_list.txt", "w")
alignlist = list(align)
cdrlist = list(cdr)
fastalist = list(fasta)
align.close()
cdr.close()
fasta.close()

ids=[]
for item in fastalist:
	if (item[0] == ">"):
		s = item.split(" ")[0]
		ids.append(s[1:len(s)])

number = len(alignlist)
for i in range(1, number):
	result.write(str(i-1)+"/"+ids[i-1]+"/")
	writeIn(result, processLine(alignlist[i]))
	result.write(cdrlist[i].split("\t")[1])

result.close()




