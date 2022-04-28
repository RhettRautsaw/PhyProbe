#!/usr/bin/env python3

import sys, os, shutil, errno, itertools
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import multiprocessing as mp
import glob
from tqdm import tqdm
import pandas as pd
try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
except:
	print("Error: biopython module is not properly installed.")
	quit()

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

CombineSources.py

Did you pull BUSCOs, Sequence Capture Loci (AHE, UCE), and other loci and are afraid there are duplicates sequences across the different datasets?
Well that's a well-placed fear since things like UCEs can often be nested within BUSCOs.

Luckily this is easy to solve with a little biopython bioinformagic. 

Before running this script you should prepare set of unaligned fastas with the format:
	>locus|sample

If your loci are not in this format, but you have a different delimiter that can separate locus from sample information, then you can use SeqReformat.py to help you.

Prepare at least 2 fastas of the options below:
	1. All named loci together
	3. All unnamed loci together

With these fastas, I will cluster sequences with either cd-hit or orthofinder (based on your preference) and preferentially keep the longest named fasta.
However, I will make sure to keep the names of the named loci associated with this sequence.

""")

parser.add_argument("-n","--named",
					type=str,
					default="",
					help="Comma separated list of fastas with named sequences")
parser.add_argument("-u", "--unnamed",
					type=str,
					default="",
					help="Comma separated list of fastas with unnamed sequences")
parser.add_argument("-i","--identity",
					type=float,
					default=0.99,
					help="cd-hit-est sequence identity threshold. (default: %(default)s)")
parser.add_argument("--orthofinder",
					action="store_true",
					default=False,
					help="Performed orthologous grouping with orthofinder instead of cd-hit-est")
parser.add_argument("-I","--inflation",
					type=float,
					default=1.5,
					help="Orthofinder inflation parameter. Only necessary if --orthofinder is used. (default: %(default)s)")
parser.add_argument("-c","--cpu",
					type=int,
					default=mp.cpu_count(),
					help="Number of cpus to be used. (default: All cpus = %(default)s)")
parser.add_argument("-o","--output",
					type=str,
					default="CombineSources",
					help="Folder name in which to do analyses and output results")
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

def read_clstr(cluster):
		# parse through the .clstr file and create a dictionary
		# with the sequences per cluster
		# open the cluster file and set the output dictionary
		cluster_file, cluster_dic = open(cluster), {}
		# parse through the cluster file and store the cluster name + sequences in the dictionary
		cluster_groups = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))
		for cluster in cluster_groups:
			name = cluster.__next__().strip()
			seqs = [seq.split('>')[1].split('...')[0] for seq in cluster_groups.__next__()]
			cluster_dic[name] = seqs
		# return the cluster dictionary
		return cluster_dic

########################################
################# SETUP ################
########################################

output=args.output
mkdir_p(output)

if args.named!=None:
	nfastas_files=[os.path.abspath(x) for x in args.named.split(",")]
	nfastas_names=[os.path.basename(x).split(".")[0] for x in nfastas_files]
	for i in nfastas_files:
		sp.call("cp " + i + " " + output, shell=True)

if args.unnamed!=None:
	unfastas_files=[os.path.abspath(x) for x in args.unnamed.split(",")]
	unfastas_names=[os.path.basename(x).split(".")[0] for x in unfastas_files]
	for i in unfastas_files:
		sp.call("cp " + i + " " + output, shell=True)

identity = args.identity * 100
inflation = args.inflation
cpus = args.cpu

print("""

CombineSources

""")

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting CombineSources :::")
print("\tNamed Sequences -> "+ ",".join(nfastas_names))
print("\tUnnamed Sequences -> "+ ",".join(unfastas_names))
if args.orthofinder:
	print("\tGrouping with --> orthofinder with inflation parameter " + str(inflation))
else:
	print("\tGrouping with --> cd-hit-est @ " + str(identity) + "%")
print("\tOutput -> " + args.output)
print("\tCPUs -> "+ str(cpus))

########################################
################# CODE #################
########################################

if args.orthofinder:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Running Orthofinder :::')
	sp.call('orthofinder -f ' + output + ' -d -og -o tmp_orthofinder -t ' + str(cpus) + ' -I ' + str(inflation), shell=True)
	sp.call('mv tmp_orthofinder/Results* ' + output + '/OF_Results', shell=True)
	sp.call('rm -r tmp_orthofinder', shell=True)
	
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Parsing Orthofinder Results :::')
	os.chdir(output)
	sp.call('rm *.fasta *.fna *.fa', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	orthogroupSeqs=sorted(glob.glob("OF_Results/Orthogroup_Sequences/*.fa"))
	orthogroups=pd.read_csv("OF_Results/Orthogroups/Orthogroups.tsv", sep="\t")
	multicopy_orthogroups=list(orthogroups["Orthogroup"])
	
	row=0
	combined_seqs=[]
	for og in tqdm(orthogroupSeqs):
		name=os.path.basename(og).split(".fa")[0]
		tmp_seq = list(SeqIO.parse(og,'fasta'))
		if name in multicopy_orthogroups:
			df=orthogroups.loc[orthogroups['Orthogroup']==name]
			new_names=[]
			for i in nfastas_names: 								# For each column with named sequences
				if pd.notna(df[i][row]): 							# check if it is empty and if not...
					for x in df[i][row].replace(' ','').split(","): # Then split up the sequence names and for each one...
						new_names.append(x.split("|")[0]) 			# split the locus from the sample name and append it to the new_name list
		new_names=sorted(new_names)
		tmp_seq2=[]
		for seq in tmp_seq:											# For each Orthogroup, see if there is a named fasta present and throw that into a new list
			tmp=seq.name.split("|")[0]
			if tmp in new_names:
				tmp_seq2.append(seq)
		max_len=0
		if tmp_seq2==[]:											# If there is no named sequences, then pick the longest unnamed sequence as the final sequence
			for seq in tmp_seq:
				if len(seq) > max_len:
					max_len = len(seq)
					final_seq = seq
					#sample=final_seq.name.split("|")[-1]
					#final_seq.name=final_seq.id=final_seq.description=name+"|"+sample
		else:
			for seq in tmp_seq2:									# If there is a named sequence, then pick the longest named sequence as the final sequence and append orthologous named sequences
				if len(seq) > max_len:
					max_len = len(seq)
					final_seq = seq
					sample=final_seq.name.split("|")[-1]
					final_seq.name=final_seq.id=final_seq.description=';'.join(new_names)+"|"+sample
		
		combined_seqs.append(final_seq)
		row+=1
	
	SeqIO.write(combined_seqs, "CombinedSources.fasta", "fasta")

##################################################################################################################
##################################################################################################################
##################################################################################################################

else:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Running cd-hit-est :::')
	os.chdir(output)
	sp.call('cat *.fasta > tmp.fasta', shell=True)
	sp.call('cd-hit-est -i tmp.fasta -o tmp_clust.fasta -c ' + str(args.identity) + ' -d 0 -sf 1 -sc 1 -T ' + str(cpus), shell=True)
	
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Parsing cd-hit Results :::')
	cluster_file=read_clstr("tmp_clust.fasta.clstr")
	cluster_fasta=list(SeqIO.parse("tmp_clust.fasta",'fasta'))
	new_fasta=[]
	for cluster in tqdm(cluster_file):
		cluster_seqs=cluster_file[cluster]
		if len(cluster_seqs)>1:
			new_names=[]
			for seq in cluster_seqs:
				new_names.append(seq.split("|")[0])
			for seq in cluster_fasta:
				if seq.id in cluster_seqs:
					sample=seq.name.split("|")[-1]
					seq.name=seq.id=seq.description=';'.join(sorted(new_names))+"|"+sample
					new_fasta.append(seq)
		else:
			new_name=''.join(cluster_seqs)
			for seq in cluster_fasta:
				if seq.id == new_name:
					new_fasta.append(seq)
	
	sp.call("rm *.fasta", shell=True)
	SeqIO.write(new_fasta, "CombinedSources.fasta", "fasta")
