#! /usr/bin/env python
# -*- coding:UTF-8 -*-

"""
------------------------------------------------------------------------------------------
It is used to get the longest mRNA of gene and filter pep and cds result with a certain length; 
usage: python xx.py -i GFF -f FASTA -g out.gff -o out.fa -s species_short_name -l 30

newly revised by Frank in Feb 22th 2020.
newly revised by Frank in Aug 15th 2020, 
	a. solve phase coding problem; 
	b. solve translate table problem.
	c. solve columns error of gff.
------------------------------------------------------------------------------------------
"""

__author__ = "Frank"
__version__ = "0.4.3"

import sys
import os
import gzip
import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from collections import OrderedDict
from optparse import OptionParser

def readfile(filein,mode):
	try:
		if filein.endswith('gz'):
			fx = gzip.open(filein,mode+'b')
		else:
			fx = file(filein,mode)
	except:
		sys.stderr.write('ERROR: fail to the IO file: %s!\n'%filein)
	
	return fx

def calculateDist(listx):
	dists = []
	for i in xrange(len(listx)):
		pos1 = int(listx[i].split('_')[0])
		pos2 = int(listx[i].split('_')[1])
		dist = pos2-pos1+1
		dists.append(dist)
	sumx = float(np.sum(dists))	
	
	return sumx

def gffAttr(idx,key_len=3,sep=';'):
	if sep in idx:
		out = idx.split(sep)[0][key_len:]
	else:
		out = idx[key_len:]
	
	return out

def reParentAttr(line):
	if 'Parent=' in line:
		parent = re.search(r'Parent=\S+;?',line).group()
		par_info = gffAttr(parent,key_len=7)
	else:
		par_info = reIDAttr(line)
	
	return par_info

def reIDAttr(line):
	idx = re.search(r'ID=\S+;?',line).group()
	id_info = gffAttr(idx)

	return id_info

def reIDAttr2(line):
	idx = re.search(r'ID=\S+;?',line)
	if idx == None:
		print line
		id_info = None
	else:
		idx = idx.group()
		id_info = gffAttr(idx)
	
	return id_info

def reNameAttr(line):
	name = re.search(r'Name=\S+;?',line).group()
	name_info = gffAttr(name,key_len=5)
	
	return name_info

def outgffAttr(info,line,split_sym=';'):
	linetype = info.split('/')[1]
	if linetype == 'gene':
		geneid = reIDAttr(line)
		geneid2 = outgene(geneid)
		line = '\t'.join(line.strip().split('\t')[0:-1])+'\tID='+geneid2+';\n'
	#elif linetype == 'mRNA' or linetype == 'CDS':
	else:
		mRNAparid = reParentAttr(line)
		mRNAparid2 = outgene(mRNAparid)
		if 'ID=' in line:
			mRNAid = reIDAttr(line)
			mRNAid2 = outgene(mRNAid)
			line = '\t'.join(line.strip().split('\t')[0:-1])+'\tID='+mRNAid2+';Parent='+mRNAparid2+';\n'
		else:
			line = '\t'.join(line.strip().split('\t')[0:-1])+'\tID='+mRNAparid2+';Parent='+mRNAparid2+';\n'

	return line

def sortedDict(hgff):
	for gene in hgff:
		if len(hgff[gene]) == 0:
			del hgff[gene]
			continue
		distcal = {}
		for x in hgff[gene]:
			distx = calculateDist(hgff[gene][x])
			distcal[x] = distx
		sorted_dist = sorted(distcal.items(), key = lambda d:d[1],reverse=True)
		sorted_key = sorted_dist[0][0]
		strand = gene.split('/')[-1]
		if strand == '+':
			sorted_value = sorted(hgff[gene][sorted_key],key=lambda x:int(x.split('_')[0]))
		else:
			sorted_value = sorted(hgff[gene][sorted_key],key=lambda x:int(x.split('_')[0]),reverse=True)
		hgff[gene] = {sorted_key:sorted_value}
	
	return hgff

def outgff(hgff,hgeneinfo,hinfo,split_sym='_'):
	houtinfo = OrderedDict()
	for gene in hgff:
		if gene in hgeneinfo:
			chrx = hgeneinfo[gene].split('/')[0]
			geneinfox = chrx+'/'+'gene'+'/'+gene
			houtinfo[geneinfox] = hinfo[geneinfox]
			for x in hgff[gene]:
				rnainfox = chrx+'/'+'mRNA'+'/'+x
				houtinfo[rnainfox] = hinfo[rnainfox]
				for i in xrange(len(hgff[gene][x])):
					cdsinfox = chrx+'/CDS/'+split_sym.join(hgff[gene][x][i].split(split_sym)[0:-1])+'/'+x
					houtinfo[cdsinfox] = hinfo[cdsinfox]
		else:
			continue

	return houtinfo

def outseq_trans_strand(seq,strand_dir):
	if strand_dir == '-':
		seqout = ''.join(list(Seq(seq, IUPAC.unambiguous_dna).reverse_complement())).upper()
	else:
		seqout = seq.upper()
	
	return seqout

def outgene(gene):  # update 2020.01.08
	#if ':' in gene:
	#	geneout = gene.split(':')[-1]
	#else:
	#	geneout = gene
	gene = gene.split(':')[-1]
	geneout = re.sub('[^.\w+]','_',gene)

	return geneout


def longestmRNA(fgff,gene_region='CDS'):
	hgff = OrderedDict()
	hmRNA_gene = {}
	hgeneinfo = {}		#{geneid1:'chr1/+',geneid2:'chr1/-'}
	hinfo = {}			#{'chr1/gene/geneid':line,'chr1/mRNA/rna_id':line'}

	for line in fgff:
		if line.startswith('#'):continue
		if len(line.strip()) == 0:continue
		gffline = line.strip().split('\t')
		geneinfo = gffline[0]+'/'+gffline[6]
		
		if gffline[2] == 'gene':				#确保gene行都有biotype字符串或都没有才行
			if 'biotype=protein_coding' in gffline[8]:
				geneid = reIDAttr(gffline[8])+'/'+gffline[6]
				hgff.setdefault(geneid,{})
				hgeneinfo[geneid] = geneinfo
				geneinfox = gffline[0]+'/'+'gene'+'/'+geneid
				hinfo[geneinfox] = line
			elif 'biotype=' not in gffline[8]:
				geneid = reIDAttr(gffline[8])+'/'+gffline[6]
				hgff.setdefault(geneid,{})
				hgeneinfo[geneid] = geneinfo
				geneinfox = gffline[0]+'/'+'gene'+'/'+geneid
				hinfo[geneinfox] = line
			else:
				continue
		else:
			if gffline[2] == 'mRNA' or gffline[2] == 'transcript':
				rna_parent = reParentAttr(gffline[8])+'/'+gffline[6]
				rna_id = reIDAttr(gffline[8])
				hmRNA_gene[rna_id] = rna_parent
				if rna_parent in hgff:
					hgff[rna_parent][rna_id] = []
				else:
					continue
				rnainfox = gffline[0]+'/'+'mRNA'+'/'+rna_id
				hinfo[rnainfox] = line
			elif gffline[2] == gene_region:
				cds_parent = reParentAttr(gffline[8])
				posinfo = gffline[3]+'_'+gffline[4]+'_'+gffline[7]
				if cds_parent in hmRNA_gene:
					genex = hmRNA_gene[cds_parent]
					if genex in hgff:
						hgff[genex][cds_parent].append(posinfo)
					else:
						continue
				else:
					continue
				cdsinfox = gffline[0]+'/CDS/'+gffline[3]+'_'+gffline[4]+'/'+cds_parent
				hinfo[cdsinfox] = line
			else:
				continue

	hgff = sortedDict(hgff)
	houtinfo = outgff(hgff,hgeneinfo,hinfo,split_sym='_')

	del hinfo
	return hgff,hgeneinfo,houtinfo

#def outfa_parse(hout,hgeneinfo,fasta,fa_out,sp_abbr="None",outseqname_type="gene"):
def outfa_parse(hout,hgeneinfo,fasta,fa_out,fp_out,sp_abbr="None",outseqname_type="gene",table=1,cds_fetch_position='True'):
	for seq in SeqIO.parse(fasta,"fasta"):
		for gene in hout:
			if len(hout[gene]) == 0:continue
			if gene in hgeneinfo:
				chrom = hgeneinfo[gene].split('/')[0]
				sym = hgeneinfo[gene].split('/')[-1]
				seqouts = []
				if chrom == seq.id:
					for x in hout[gene]:
						phase_list = []
						for i in xrange(len(hout[gene][x])):
							pos1 = int(hout[gene][x][i].split('_')[0])
							pos2 = int(hout[gene][x][i].split('_')[1])
							seqout = str(seq.seq)[pos1-1:pos2]
							seqouts.append(seqout)
							phase_list.append(hout[gene][x][i].split('_')[-1])
						
					del hout[gene]
					#del hgeneinfo[gene]
					if sym == '+':
						outseqs = ''.join(seqouts)
					else:
						outseqs = ''.join(seqouts[::-1])
					seq_ret = outseq_trans_strand(outseqs,sym)
					
					phase = phase_list[0]
					if phase == '.':
						seq_retx = seq_ret
					else:
						seq_retx = seq_ret[int(phase):]
					
					gene = outgene(gene.split('/')[0])
					x = outgene(x)
					if outseqname_type == "gene":
						if sp_abbr == "None":
							if cds_fetch_position == 'True':
								fa_out.write('>'+gene+'\n'+seq_ret+'\n')
							else:
								fa_out.write('>'+gene+'\n'+seq_retx+'\n')
							fp_out.write('>'+gene+'\n'+str(Seq(seq_retx).translate(table=table))+'\n')
						else:
							if cds_fetch_position == 'True':
								fa_out.write('>'+gene+'_'+sp_abbr+'\n'+seq_ret+'\n')
							else:
								fa_out.write('>'+gene+'_'+sp_abbr+'\n'+seq_retx+'\n')
							fp_out.write('>'+gene+'_'+sp_abbr+'\n'+str(Seq(seq_retx).translate(table=table))+'\n')
					elif outseqname_type == "mRNA":
						if sp_abbr == "None":
							if cds_fetch_position == 'True':
								fa_out.write('>'+x+'\n'+seq_ret+'\n')
							else:
								fa_out.write('>'+x+'\n'+seq_retx+'\n')
							fp_out.write('>'+x+'\n'+str(Seq(seq_retx).translate(table=table))+'\n')
						else:
							if cds_fetch_position == 'True':
								fa_out.write('>'+x+'_'+sp_abbr+'\n'+seq_ret+'\n')
							else:
								fa_out.write('>'+x+'_'+sp_abbr+'\n'+seq_retx+'\n')
							fp_out.write('>'+x+'_'+sp_abbr+'\n'+str(Seq(seq_retx).translate(table=table))+'\n')
					else:
						sys.stderr.write('[ERROR] Something wrong with the -t option, please check it.')
				else:
					continue
			else:
				continue

def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option('-i',"--input-gff",action="store",type="string",dest="input_gff",help="Annotation file in gff format. [required]")
	parser.add_option('-f',"--input-fasta",action="store",type="string",dest="input_fasta",help="genome file in fasta format. [required]")
	parser.add_option('-g',"--out-gff",action="store",type="string",dest="output_gff",help="output files(s) in gff format with longest mRNA. [required]")
	parser.add_option('-o',"--out-fa-prefix",action="store",type="string",dest="output_fasta",help="prefix of output files(s) in fasta format. [required]")
	parser.add_option('-s',"--species-shortname",action="store",type="string",dest="species_shortname",help="add species shortname in sequence name",default="None")
	parser.add_option('-r',"--region",action="store",dest="gene_region",help="the region of genes, you can choose 'CDS' or 'exon'",type="choice",choices=["CDS","exon"],default="CDS")
	parser.add_option('-t',"--type",action="store",dest="outseqname_type",help="type of out sequence name,you can choose 'gene' or 'mRNA'",type="choice",choices=["gene","mRNA"],default="mRNA")
	parser.add_option('-l',"--length",action="store",dest="pep_length",help="the filtered length of protein sequence",type="int",default=30)
	parser.add_option('--tt',action='store',dest="transl_table",help="codon table used for translating CDS sequence, you can learn more from the link: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi",type="int",default=1)
	parser.add_option('--cfp',action='store',dest='cds_fetch_position',help='CDS sequences extracted by position or not, True or False, default="False"',type='choice',choices=["True","False"],default="False")

	(options,args)=parser.parse_args()

	if not (options.input_fasta and options.input_gff):
		parser.print_help()
		sys.exit(0)
	
	with readfile(options.input_gff,'r') as fgff:
		hout,hgeneinfo,houtinfo = longestmRNA(fgff,gene_region=options.gene_region)
	with readfile(options.output_gff,'w') as fgff_out:
		for outinfo in houtinfo:
			fgff_out.write(outgffAttr(outinfo,houtinfo[outinfo]))

	with readfile(options.input_fasta,'r') as fasta, readfile(options.output_fasta+'.cds','w') as fa_out, readfile(options.output_fasta+'.pep','w') as fp_out:
		outfa_parse(hout,hgeneinfo,fasta,fa_out,fp_out,sp_abbr=options.species_shortname,outseqname_type=options.outseqname_type,table=options.transl_table,cds_fetch_position=options.cds_fetch_position)
	os.system('python pep_cds_filter.py %s %s %s'%(os.path.abspath(options.output_fasta+'.pep'),os.path.abspath(options.output_fasta+'.cds'),str(options.pep_length)))

if __name__ == "__main__":
	main()
