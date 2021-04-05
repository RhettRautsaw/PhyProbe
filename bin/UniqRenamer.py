#!/usr/bin/env python3

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import gzip
import copy
from tqdm import tqdm
from p_tqdm import p_map
from dfply import *
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
	import glob
except:
	print("Error: glob is not properly installed.")
	quit()

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

Check sequences against database and rename appropriately.

Assumes targets have the format:
>target_name|sample_name

But delimiter | can be changed with the -d flag.

""")

########################################
############### ARGUMENTS ##############
########################################

parser.add_argument("-i","--input",
					type=str,
					help="File with sequences to check agains target database (default: %(default)s)")
parser.add_argument("-t","--targets",
					type=str,
					help="Target fasta to annotate orthogroups (to ensure proper duplicate removal) (default: %(default)s)")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Delimiter used to separate target/gene names from sample names (default: %(default)s)")
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

########################################
################# SETUP ################
########################################

input=os.path.abspath(args.input)
folder=input.split("/")[:-1]
folder='/'.join(folder)
targets_name = os.path.abspath(args.targets)
delim = args.delim
cpus = args.cpu

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting UniqRenamer...")
print("\tInput -> "+ input)
print("\tTargets Filter -> "+ targets_name)
print("\tDelim -> "+ delim)
print("\tThreads -> "+ str(cpus))

########################################
################# CODE #################
########################################

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Reading ' + args.input)
sequences = list(SeqIO.parse(input,"fasta"))

if not os.path.isfile(targets_name + ".nin"):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running makeblastdb :::\n")
	sp.call('makeblastdb -in ' + targets_name + ' -dbtype nucl', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

if not os.path.isfile(folder + '/blast_results.tsv'):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running BLAST :::\n")
	sp.call('blastn -query ' + input + ' -db ' + targets_name + ' -outfmt 6 -num_threads ' + str(cpus) + ' -evalue 0.0001 -out ' + folder + '/blast_results.tsv', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL) #-max_target_seqs 1 -max_hsps 1

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Renaming sequences based on BLAST results :::\n")
blastfile = pd.read_csv(folder + '/blast_results.tsv', sep="\t", names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
blastfile = blastfile >> mutate(tmp = X.sseqid) >> separate(X.tmp, ['ref'], sep="\\"+delim) >> mutate(tmp = X.qseqid) >> separate(X.tmp, ['qseqid2'], sep="\\"+delim)

#blastfile = blastfile[['qseqid2','ref']].drop_duplicates()
#print(blastfile.to_string(index=False, header=None)+'\n')
#
#dct = {}
#for i, j in zip(list(blastfile.qseqid2), list(blastfile.ref)):
#	dct.setdefault(i, []).append(j)
#
#sequences2=[]
#for seq in tqdm(sequences):
#	locus=seq.name.split(delim)[0]
#	sample=seq.name.split(delim)[1]
#	if locus in dct:
#		for ref in dct[locus]:
#			seq2 = copy.deepcopy(seq)
#			seq2.name = seq2.id = seq2.description = ref + delim + sample + delim + locus
#			sequences2.append(seq2)
#	else:
#		sequences2.append(seq)
#

blastfile2 = blastfile >> group_by(X.qseqid, X.ref) >> summarize(qstart_min=colmin(X.qstart), qend_max=colmax(X.qend))
match_list = list(blastfile2.qseqid)

sequences2=[]
for seq in tqdm(sequences):
	if seq.name in match_list:
		locus=seq.name.split(delim)[0]
		sample=seq.name.split(delim)[1]
		tmp_df = blastfile2.loc[blastfile2['qseqid'] == seq.name]
		for i in range(len(tmp_df)):
			ref = ''.join(list(tmp_df.iloc[[i]].ref))
			start = int(tmp_df.iloc[[i]].qstart_min)-1
			end = int(tmp_df.iloc[[i]].qend_max)-1
			seq2 = copy.deepcopy(seq)
			seq2.name = seq2.id = seq2.description = ref + delim + sample + delim + locus
			seq2.seq = seq.seq[start:end]
			sequences2.append(seq2)
	else:
		sequences2.append(seq)

sequences = sequences2

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Writing results :::")
handle=open(input.split('.')[0]+"_renamed.fasta","w")
writer = FastaIO.FastaWriter(handle)
writer.write_file(sequences)
handle.close()


