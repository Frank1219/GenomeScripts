#! /usr/bin/env python
# -*- coding:UTF-8 -*-

__author__ = 'Frank' 

import sys
import os
import re
from fa_cds_longest import readfile
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC

def filter_IDname(id):
	if ':' in id and id.split(':')[0] == 'gene':
		id2 = id.strip().split(':')[1]
	elif '_evm' in id:
		id2 = id.strip().split('_evm')[0]+'.0'+'_'+id.strip().split('_')[-1]
	else:
		id2 = id
	
	return id2

def Seqfilter(fpep,fcds,filter_leng,fout_pep,fout_cds):
	hseq = {}
	for seq in SeqIO.parse(fpep,"fasta"):
		id = seq.id.split()[0]
		#id2 = filter_IDname(id)
		if ('X' not in str(seq.seq)) and ('*' not in str(seq.seq)[:-1]) and len(str(seq.seq)) >= filter_leng:
			#fout_pep.write('>'+id2+'\n'+str(seq.seq)+'\n')
			fout_pep.write('>'+id+'\n'+str(seq.seq)+'\n')
			hseq[id] = None
		else:
			continue

	for seq in SeqIO.parse(fcds,"fasta"):
		#idx = filter_IDname(seq.id)
		#if idx in hseq:
		if seq.id in hseq:
			fout_cds.write('>'+seq.id+'\n'+str(seq.seq)+'\n')
		else:
			continue

def main():
	filter_leng = int(sys.argv[3])
	with readfile(sys.argv[1],'r') as fpep, readfile(sys.argv[2],'r') as fcds, readfile(os.path.basename(sys.argv[1])+'.filter','w') as fout_pep, readfile(os.path.basename(sys.argv[2])+'.filter','w') as fout_cds:
		Seqfilter(fpep,fcds,filter_leng,fout_pep,fout_cds)

if __name__ == "__main__":
	main()
