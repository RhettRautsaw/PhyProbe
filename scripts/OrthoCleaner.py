#!/usr/bin/env python

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import gzip
import multiprocessing as mp
from p_tqdm import p_map
import pandas as pd
import numpy as np
import glob

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

Extract single copy loci from OrthoFinder (including loci with missing individuals).

Results will have the format:
>locus_name|sample_name

But delimiter | can be changed with the -d flag.

""")

########################################
############### ARGUMENTS ##############
########################################

parser.add_argument("-f","--folder",
					type=str,
					default='00_ortho_res',
					help="Base folder with orthofinder results (default: %(default)s)")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Delimiter used to separate target/gene names from sample names (default: %(default)s)")
parser.add_argument("-p","--percent",
					type=int,
					default=50,
					help="Percent missing data allowed. (default: %(default)s)")
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

def OG_extract(file):
	file_name = file.split("/")[-1]
	prefix = file_name.split(".fa")[0]
	if prefix in Orthogroups:
		sp.call('cp '+ file + ' ' + folder + '/OrthoCleaner/' + prefix + '.fasta',shell=True)
		sp.call("perl -pi -e 's/>/>" + prefix + "\\" + delim + "/g' " + folder + "/OrthoCleaner/" + prefix + '.fasta', shell=True)
		#sp.call("perl -pi -e 's/>/>" + prefix + "\\"+ delim + "/g' "+ folder + "/OrthoCleaner/" + prefix + '.fasta', shell=True)
		

########################################
################# SETUP ################
########################################

if args.folder==None:
	print("No OrthoFinder output folder specified. Use -f/--folder flag to specify input")
	quit()
else:
	folder = os.path.abspath(args.folder)

delim = args.delim
percent=args.percent/100
cpus = args.cpu

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting OrthoCleaner...")
print("\tInput:\t\t"+ folder)
print("\tDelim:\t\t"+ delim)
print("\tMissing:\t<"+ str(args.percent) + "%")
print("\tCPUs:\t\t"+ str(cpus))

########################################
################# CODE #################
########################################

file = glob.glob(folder + '/Results_*/Orthogroups/Orthogroups.GeneCount.tsv')[0]
GeneCount = pd.read_csv(file, sep='\t', index_col=0)
total_samps = GeneCount.shape[1]-1

# Determine if the max gene count is less than or equal to 1 (single copy genes only)
GeneCount2 = GeneCount.drop(columns=['Total'])
GeneCount['sub1'] = list(GeneCount2.max(axis=1) <= 1)

# Determine if the total gene count has less than X missing data
per_miss = round(percent * total_samps)
GeneCount['sub2'] = list(GeneCount.Total >= per_miss)

# Subset Orthogroups
GeneCount = GeneCount[GeneCount["sub1"]]
GeneCount = GeneCount[GeneCount["sub2"]]

Orthogroups = list(GeneCount.index.values)

mkdir_p(folder + '/OrthoCleaner')

files=glob.glob(folder + '/Results_*/Orthogroup_Sequences/*.fa')
print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Extracting ' + str(len(Orthogroups)) + ' single-copy orthogroups of ' + str(len(files)) + ' total orthogroups :::')

results = p_map(OG_extract, files, num_cpus=cpus)

sp.call("cat " + folder + "/OrthoCleaner/*.fasta > " + folder + "/OrthoCleaner.fasta", shell=True)
sp.call("perl -pe 's/\\" + delim + ".*//g' " + folder + "/OrthoCleaner.fasta > " + folder + "/OrthoCleaner_renamed.fasta", shell=True)