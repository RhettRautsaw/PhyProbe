#!/usr/bin/env python3

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import glob
from dfply import *
from p_tqdm import p_map
from collections import defaultdict
import pyfastx

try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.SeqIO import FastaIO
except:
	print("Error: biopython module is not properly installed.")
	quit()

try:
	import numpy as np
except:
	print("Error: numpy is not properly installed.")
	quit()

try:
	import pandas as pd
except:
	print("Error: pandas is not properly installed.")
	quit()

try:
	from dfply import *
except:
	print("Error: dfply is not properly installed.")
	quit()

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""s

LocusExtractor
1. Parse reads to target loci with blast or minimap2
2. Assemble contigs with spades
3. Extract consensus sequence (sam2con.py)

Results will have the format:
>locus_name|sample_name

Delimiter | can be changed with the -d flag.

""")

parser.add_argument("-r","--reads",
					nargs='+',
					help="Space-separated list of reads in fastq(.gz) format (default: %(default)s)")
parser.add_argument("-t","--targets",
					type=str,
					help="Reference fasta with target loci (default: %(default)s)")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Delimiter used to separate target/gene names from sample names (default: %(default)s)")
parser.add_argument("-p","--prefix",
					type=str,
					help="Unique name of sample for output (default: %(default)s)")
parser.add_argument("-o","--output",
					type=str,
					default="06_targets",
					help="Folder in which to export results (default: %(default)s)")
parser.add_argument("--percent",
					type=int,
					default=25,
					help="Percent missing data (gaps) allowed in a locus. (default: %(default)s)")
parser.add_argument("--blast",
					action="store_true",
					default=False,
					help="Use blast instead of minimap2 mapping")
parser.add_argument("-c","--cpu",
					type=int,
					default=8,
					help="Number of threads to be used in each step. (default: %(default)s)")
args=parser.parse_args()

########################################
############### FUNCTIONS ##############
########################################

def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc: # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else: raise

def Parser(ref):
	folder = os.path.join(genes_folder,ref)
	file = os.path.join(folder,ref)
	mkdir_p(folder)
	
	reads_list = pd.DataFrame(read_dict[ref])
	reads_list.to_csv(file+"_reads.list",header=False, index=False)
	
	text_file = open(file+"_ref.list", "w")
	n = text_file.write(best_refs[ref])
	text_file.close()
	
	sp.call("seqtk subseq " + sub_reads + " " + file + "_reads.list > " + file + "_reads.fasta", shell=True)
	sp.call("seqtk subseq " + targets_name + " " + file + "_ref.list > " + file + "_ref.fasta", shell=True)

def Assembler(ref):
	folder = os.path.join(genes_folder,ref) # output/genes/ref
	file = os.path.join(folder, ref) # output/genes/ref/ref
	
	sp.call("spades.py --only-assembler --threads 1 --cov-cutoff 8 -s " + file + "_reads.fasta -o " + os.path.join(folder,"spades"), shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	if os.path.isfile(os.path.join(folder,"spades/contigs.fasta")):
		sp.call("mv " + os.path.join(folder,"spades/contigs.fasta") + " " + file + "_contigs.fasta", shell=True)
	else:
		sp.call("mv " + file + "_reads.fasta" + " " + file + "_contigs.fasta", shell=True)
	sp.call('minimap2 -ax splice ' + file + "_ref.fasta " + file + "_contigs.fasta > " + os.path.join(folder,"align.sam"), shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call('sam2con.py -i ' + os.path.join(folder,"align.sam") + ' -o ' + folder + ' -p "' + ref + '" -p2 "' + ref + '|' + prefix + '"', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	if os.path.isfile(file + '_consensus.fasta'):
		tmp_seq = list(SeqIO.parse(file + '_consensus.fasta','fasta'))
		for seq in tmp_seq:
			c = seq.seq.count("-")
			l = len(seq.seq)
			if c/l < percent:
				output_handle = open(folder + ".fasta",'a')
				SeqIO.write(seq, output_handle, "fasta")
				output_handle.close()
				sp.call('cat '+ folder + '.fasta >> ' + os.path.join(output,prefix) + '.targets.fasta', shell=True)
			else:
				sp.call('rm -r ' + os.path.join(folder), shell=True)
	else:
		sp.call('rm -r ' + os.path.join(folder), shell=True)


########################################
################# SETUP ################
########################################

reads = [os.path.abspath(x) for x in args.reads]
targets_name = os.path.abspath(args.targets)
delim = args.delim
prefix = args.prefix
output = args.output
percent = args.percent/100
blast = args.blast
cpus = args.cpu

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting LocusExtractor...")
start_dir = os.getcwd()
print("\tReads -> "+ ' '.join(reads))
print("\tTargets -> "+ targets_name)
print("\tDelimiter -> "+ delim)
print("\tPrefix -> "+ prefix)
print("\tOutput -> "+ output)
print("\tPercent -> "+ str(percent))
print("\tThreads -> "+ str(cpus))
print("\tUse BLAST -> "+ str(blast))

########################################
################# CODE #################
########################################

prefix2 = os.path.join(output,prefix)
mkdir_p(output)
genes_folder=os.path.join(output,"genes")
mkdir_p(genes_folder)

targets = pyfastx.Fasta(targets_name)
all_refs = set([seq.name.split(delim)[0] for seq in targets])

all_reads = prefix2 + ".combined.fasta"
sub_reads=prefix2 + ".mapped.fasta"


# CONCATENATE READS TOGETHER
if os.path.isfile(all_reads):
	pass
else:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Concatenating Read Files :::')
	for file in reads:
		sp.call("seqtk seq -a " + file + " >> " + all_reads, shell=True)


# ALIGN READS TO TARGET FILE
if blast:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Aligning Reads to Targets :::')
	sp.call("blastn -query " + all_reads + " -db " + targets_name + " -max_target_seqs 10 -evalue 1e-10 -outfmt 6 -num_threads " + str(cpus) + " > " + prefix2 + ".blast", shell=True)
	alignfile = pd.read_csv(prefix2 + ".blast", sep="\t", names=['qseqid', 'sseqid', 'length', 'score'], usecols=(0,1,3,11))
	sp.call("cut -f1 " + prefix2 + ".blast | sort | uniq > " + prefix2 + "_map.list", shell=True)
else:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Aligning Reads to Targets :::')
	sp.call("minimap2 -x sr -t " + str(cpus) + " " + targets_name + " " + all_reads + " > " + prefix2 + ".paf", shell=True)
	alignfile = pd.read_csv(prefix2 + ".paf", sep="\t", names=['qseqid', 'sseqid', 'length', 'score'], usecols=(0,5,6,11))
	sp.call("cut -f1 " + prefix2 + ".paf | sort | uniq > " + prefix2 + "_map.list", shell=True)


# SUBSAMPLE READS WITH HITS FOR SPEED
print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Subsampling Reads with Hits :::')
sp.call("seqtk subseq " + all_reads + " " + prefix2 + "_map.list > " + sub_reads, shell=True)


# COUNT HITS
print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Counting # of loci that have hits :::')
alignfile = alignfile >> mutate(tmp = X.sseqid) >> separate(X.tmp, ['ref'], sep="\\"+delim)
best_refs = alignfile >> group_by(X.ref) >> arrange(X.length, X.score, ascending=False) >> head(1)
best_refs=dict(zip(list(best_refs.ref),list(best_refs.sseqid)))
tmp = alignfile[['ref','qseqid']].drop_duplicates()
read_dict = defaultdict(list)
for k, v in zip(list(tmp.ref), list(tmp.qseqid)):
	read_dict[k].append(v)

refs = list(alignfile.ref.unique())
print("\t " + str(len(alignfile)) + " hits for " + str(len(refs)) + "/" + str(len(all_refs)) + " targets")

del alignfile
del tmp

# RUN PARSER
print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Parsing Reads and References to Folders :::')
results=p_map(Parser, refs, num_cpus=cpus)


# RUN SPADES ASSEMBLY AND SCAFFOLD (REMOVING GAPS).
print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Assembling and Extracting Consensus Sequences :::')
results=p_map(Assembler, refs, num_cpus=cpus)

sp.call("rm -r " + all_reads, shell=True)
sp.call("rm -r " + sub_reads, shell=True)
sp.call("rm -r " + prefix2 + "_map.list", shell=True)
sp.call("rm -r " + prefix2 + ".paf", shell=True)
