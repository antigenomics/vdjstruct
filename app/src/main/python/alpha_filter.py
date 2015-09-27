import sys
from Bio import SeqIO
from Bio import SeqFeature

records = list(SeqIO.parse(sys.argv[1], "genbank"))
i=0
indices_list = []
while(i < len(records)):
	desc = records[i].description
	length = len(desc)
	if ((desc.find("T-cell", 0, length-1) != -1) or (desc.find("T-lymphocyte", 0, length-1) != -1) or (desc.find("TCR", 0, length-1) != -1)):
		if(records[i].annotations["references"][0].pubmed_id.isalnum()):
			if (desc.find("alpha", 0, length-1) != -1): 
				i+=1
				continue
	records.pop(i)
SeqIO.write(records, "generated/alpha_filtered.fasta", "fasta")
