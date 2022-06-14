#!/usr/bin/env python3

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import multiprocessing as mp
import glob
import pandas as pd
import numpy as np
from tqdm import tqdm
from p_tqdm import p_map
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

Blah blah blah

""")

parser.add_argument("-f","--folder",
					type=str,
					default="05_trimal",
					help="Folder of trimmed alignments")
parser.add_argument("-o","--output",
					type=str,
					default="DataSubsets",
					help="Folder for data subsets")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Delimiter used to separate target/gene names from sample names (default: %(default)s)")
parser.add_argument("-g","--gaps",
					type=int,
					default=30,
					help="Percent of sequence allowed to be gaps (only one value allowed here. Recommend to leave as default: %(default)s)")
parser.add_argument("-p1","--percent1",
					type=str,
					default="5,25,40,50,60,75,95",
					help="Percent of missing taxa allowed in a given alignment (comma separate for multiple thresholds) (default: %(default)s)")
parser.add_argument("-p2","--percent2",
					type=str,
					default="5,25,40,50,60,75,95",
					help="Percent of missing genes allowed for a given taxa (comma separate for multiple thresholds) (default: %(default)s)")
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

def MissingData(file):
	locus = file.split(".")[0].split("/")[-1]
	tmp_seq = list(SeqIO.parse(file,'fasta'))
	clean_seqs = []
	data = dict()
	for seq in tmp_seq:
		sample = seq.name.split("|")[0]
		c = seq.seq.count("-")
		l = len(seq.seq)
		completeness = 1-(c/l)
		if completeness > (1-gaps):
			clean_seqs.append(seq)
			data[sample] = completeness
	
	df = pd.DataFrame({'sample':list(data.keys()), str(locus) : list(data.values())})
	df.to_csv(os.path.join(tmpdir, locus) + ".csv", index=None)
	
	handle=open(os.path.join(tmpdir, locus) + ".fasta",'w')
	writer = FastaIO.FastaWriter(handle)
	writer.write_file(clean_seqs)
	handle.close()

def Prepare(file):
	locus = os.path.split(file)[-1].split(".")[0]
	clean_seqs = []
	if locus not in bad_loci:
		tmp_seq = list(SeqIO.parse(file,'fasta'))
		for seq in tmp_seq:
			sample = seq.name.split("|")[0]
			if sample not in bad_samples2:
				clean_seqs.append(seq)
	
	if len(clean_seqs) > 0:
		handle=open(os.path.join(resdir, locus) + ".fasta",'w')
		writer = FastaIO.FastaWriter(handle)
		writer.write_file(clean_seqs)
		handle.close()


########################################
################# SETUP ################
########################################

if args.folder==None:
	print("No input folder specified. Use -f/--folder flag to specify")
	quit()
else:
	folder = os.path.abspath(args.folder)
	files = glob.glob(folder + "/*.fasta")

if args.output==None:
	print("No output location specified. Use -o/--output flag to specify")
	quit()
else:
	output = args.output
	tmpdir = output+"/drop_gappy_samples"
	mkdir_p(output)
	mkdir_p(tmpdir)

delim=args.delim
gaps = args.gaps/100
percent1 = [ (int(x)/100) for x in args.percent1.split(",") ]
percent2 = [ (int(x)/100) for x in args.percent2.split(",") ]
cpus = args.cpu

print("""

Blah blah blah

""")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting...")
print("\tInput folder -> "+ args.folder)
print("\tOutput folder -> "+ output)
print("\tDelimiter -> " + delim)
print("\tPercent Gaps Allowed in Sequence -> " + str(args.gaps) + "%")
print("\tPercent Missing Taxa Allowed in Alignment -> " + str(args.percent1) + "%")
print("\tPercent Missing Genes Allowed in Taxa -> " + str(args.percent2) + "%")
print("\tNumber of CPU -> "+ str(cpus))

########################################
################# CODE #################
########################################

results = pd.DataFrame()

p_map(MissingData, files, num_cpus=cpus)

files2=glob.glob(tmpdir+"/*.csv")
files3=glob.glob(tmpdir+"/*.fasta")

for file in tqdm(files2):
	tmp=pd.read_csv(file)
	tmp.set_index('sample', inplace=True)
	tmp=tmp.transpose()
	results=results.append(tmp, sort=True)

results.to_csv(output+"/results.csv")

# Dropping samples with only a handful of genes
#sample_sums=results.count(axis=0)
#bad_samples=list(sample_sums[sample_sums<sample_sums.quantile(0.1)].index)
#results2=results.drop(bad_samples, axis=1)
#print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Dropping " + str(len(bad_samples)) + " outlier samples :::")
bad_samples=[]
results2=results

# Removing loci with 
locus_sums=results2.count(axis=1)
for x in percent1:
	for y in percent2:
		resdir=os.path.join(output,"Genes"+str(round(x*100))+"Taxa_Taxa"+str(round(y*100))+"Genes")
		bad_loci=list(locus_sums[(locus_sums/results2.shape[1])<x].index)
		results3=results2.drop(bad_loci, axis=0)
		
		sample_sums2=results3.count(axis=0)
		bad_samples2=list(sample_sums2[(sample_sums2/results3.shape[0])<y].index)
		results4=results3.drop(bad_samples2, axis=1)
		results4=results4.dropna(axis=0,how='all').dropna(axis=1,how='all')
		
		bad_samples2=bad_samples+bad_samples2
		
		if 0 not in results4.shape:
			print(resdir + str(results4.shape))
			
			mkdir_p(resdir)
			
			results4.to_csv(resdir+"/results.csv")
			df = pd.DataFrame(bad_loci)
			df.to_csv(resdir+"/badloci.list")
			df = pd.DataFrame(bad_samples2)
			df.to_csv(resdir+"/badsamples.list")
			
			p_map(Prepare, files3, num_cpus=cpus)
