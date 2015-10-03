#!/bin/bash

python src/main/python/beta_filter.py data/tcr.gbk
mixcr align generated/beta_filtered.fasta generated/beta_alignments.vdjca
mixcr exportAlignments generated/beta_alignments.vdjca generated/beta_alignments.txt
mixcr exportAlignments -readid -aaFeature CDR3 generated/beta_alignments.vdjca generated/beta_cdr3.txt
python src/main/python/beta_list.py

python src/main/python/alpha_filter.py data/tcr.gbk
mixcr align generated/alpha_filtered.fasta generated/alpha_alignments.vdjca
mixcr exportAlignments generated/alpha_alignments.vdjca generated/alpha_alignments.txt
mixcr exportAlignments -readid -aaFeature CDR3 generated/alpha_alignments.vdjca generated/alpha_cdr3.txt
python src/main/python/alpha_list.py
