#!/usr/bin/env python

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import multiprocessing as mp
import glob
try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.SeqIO import FastaIO
except:
	print("Error: biopython module is not properly installed.")
	quit()

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

Perform self-BLAST search and then cluster sequences with Markov Clustering (MCL)

Modes for final output: 
	--mode 1	Sequences only from cluster 1 (the largest cluster)
	--mode 2	Largest sequence from each cluster
	--mode 3	All sequences, separate fastas per cluster

""")

parser.add_argument("-i","--input",
					type=str,
					help="Input fasta")
parser.add_argument("--inflation",
					type=float,
					default=2.0,
					help="MCL Inflation Parameter (default: %(default)s)")
parser.add_argument("-m", "--mode",
					type=str,
					default="3",
					help="--mode 1	Sequences only from cluster 1 (the largest cluster)\n\
						--mode 2	Largest sequence from each cluster\n\
						--mode 3	All sequences, separate fastas per cluster\n\
						(default: %(default)s)")
parser.add_argument("-c","--cpu",
					type=int,
					default=mp.cpu_count(),
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

input = os.path.abspath(args.input)
input_name = os.path.basename(input)
input_name2 = os.path.splitext(input_name)[0]
errlog = open(input+".error", 'a')

inflation = args.inflation
mode = args.mode
cpus = args.cpu

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting MCL :::")
print("\tInput:\t\t"+ input)
print("\tMCL Inflation:\t" + str(inflation))
print("\tMode:\t\t" + mode)
print("\tCPU:\t\t"+ str(cpus))

########################################
################# CODE #################
########################################

seqs = list(SeqIO.parse(input,'fasta'))
mkdir_p("MCL_out")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running BLAST :::")
sp.call("makeblastdb -in " + input + " -dbtype nucl", shell=True, stdout=errlog, stderr=errlog)
sp.call("blastn -query " + input + " -db " + input + " -num_threads " + str(cpus) + " -evalue 1e-10 -task dc-megablast -outfmt 6 -out " + input_name + ".blast", shell=True, stdout=errlog, stderr=errlog)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running MCL :::")
sp.call("cut -f 1,2,11 " + input_name + ".blast > " + input_name + ".abc", shell=True, stdout=errlog, stderr=errlog)
sp.call("mcxload -abc " + input_name + ".abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o " + input_name + ".mci -write-tab " + input_name + ".tab", shell=True, stdout=errlog, stderr=errlog)
sp.call("mcl " + input_name + ".mci -I " + str(inflation) + " -use-tab " + input_name + ".tab", shell=True, stdout=errlog, stderr=errlog)

clusters = open(glob.glob("out." + input_name + ".mci.I*")[0], "r").read().split("\n")
clusters = [i for i in clusters if i]

if(mode == "1"):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Extracting Largest Cluster :::")
	#sp.call("head -1 out." + input_name + ".mci.I* | sed 's/\\t/\\n/g' > " + input_name + ".subset.mcl", shell=True, stdout=errlog, stderr=errlog)
	seqs2=[]
	clust1 = clusters[0].split("\t") #open(input_name + ".subset.mcl", "r").read().split("\n")
	for seq in seqs:
		if seq.name in clust1:
			seqs2.append(seq)
	
	output_handle = open("MCL_out/" + input_name2 + ".mcl.fasta", 'w')
	SeqIO.write(seqs2, output_handle, "fasta")
	output_handle.close()

if(mode=="2"):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Extracting Largest Sequence per Cluster :::")
	d={}
	for seq in seqs:
		d[seq.name]=len(seq)
	
	og = 0
	for group in clusters:
		l = 0
		tmp_clusters = group.split("\t")
		for seq in seqs:
			if seq.name in tmp_clusters and d[seq.name] > l:
				out_seq = seq
				l = len(out_seq)
		output_handle = open("MCL_out/OG" + str(og).zfill(7) + ".largest.fasta", 'w')
		SeqIO.write(out_seq, output_handle, "fasta")
		output_handle.close()
		og = og + 1

if(mode == "3"):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Separating Sequences into Clusters :::")
	og = 0
	for group in clusters:
		tmp_clusters = group.split("\t")
		out_seqs = []
		for seq in seqs:
			if seq.name in tmp_clusters:
				out_seqs.append(seq)
		
		output_handle = open("MCL_out/OG" + str(og).zfill(7) + ".fasta", 'w')
		SeqIO.write(out_seqs, output_handle, "fasta")
		output_handle.close()
		og = og + 1

sp.call("rm " + input_name + ".* out." + input_name + ".*", shell=True, stdout=errlog, stderr=errlog)