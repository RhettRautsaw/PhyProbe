#!/usr/bin/env python3

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import multiprocessing as mp
import glob
import copy
import gc
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

BuscoCleaner
1. If no fna files, extract sequences from original assembly
1. Check multi-copy busco sequences for 100% duplicates
	- Because they are 100% duplicates, they are not true multi-copy genes
	- This will result in more single-copy genes
2. Rename sequences and concatenate

Results will have the format:
>locus_name|sample_name

Delimiter | can be changed with the -d flag.

""")

parser.add_argument("-f","--folder",
					type=str,
					help="BUSCO output folder")
parser.add_argument("-n","--name",
					type=str,
					help="Name to apply to sequence")
parser.add_argument("-o","--output",
					type=str,
					help="Optional argument. For concatenating multiple busco runs, specify path here")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Delimiter used to separate locus names from sample names (default: %(default)s)")
parser.add_argument("-a","--assembly",
					type=str,
					help="Original assembly on which BUSCO was performed \n Only required if fna files were not produced (default for BUSCO v5)")
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

def Distributor(file):
	new_file=file.replace(".faa",".fna")
	new_file_seqs=[]
	locus=file.split("/")[-1].split(".")[0]
	locus_seqs=list(SeqIO.parse(file,"fasta"))
	for seq in locus_seqs:
		contig=seq.id.split(":")[0]
		start=int(seq.id.split(":")[1].split("-")[0])+1
		end=int(seq.id.split(":")[1].split("-")[1])+1
		for acontig in assembly:
			if acontig.id == contig:
				tmp=copy.deepcopy(acontig)
				tmp.seq=tmp.seq[start:end]
				tmp.id=seq.id
				new_file_seqs.append(tmp)
	SeqIO.write(new_file_seqs, new_file, "fasta")

def mcChecker(locus):
	sp.call("cd-hit-est -i " + locus + " -o " + locus + ".clust.fasta -d 0 -c 1 &> " + locus + ".clust.log", shell=True)
	sequences = list(SeqIO.parse(locus + ".clust.fasta","fasta"))
	if len(sequences)==1:
		handle=open(locus.replace("multi_copy_busco_sequences","single_copy_busco_sequences"), "w")
		writer = FastaIO.FastaWriter(handle, wrap=None)
		writer.write_file(sequences)
		handle.close()
		
		shutil.move(locus, mc_folder+"/Duplicates/")
		shutil.move(locus.replace("fna","faa"), mc_folder+"/Duplicates/")
	
	sp.call("rm " + locus + ".clust*", shell=True)

def Renamer(locus):
	locus_name=locus.split('/')[-1]
	locus_name=locus_name.split(".fna")[0]
	sp.call("sed 's/>.*/>" + locus_name + "\\" + delim + name + "/g' " + locus + " > " + os.path.join(folder,"BuscoCleaner",locus_name) + ".fasta", shell=True)
	sp.call("cat " + os.path.join(folder,"BuscoCleaner",locus_name) + ".fasta >> " + folder + "/" + name + "_busco.fasta", shell=True)

########################################
################# SETUP ################
########################################

if args.folder==None:
	print("No busco output folder specified. Use -f/--folder flag to specify input")
	quit()
else:
	folder = os.path.abspath(args.folder)

if args.name==None:
	print("No sample name specified. Use -n/--name flag to specify input")
	quit()
else:
	name = args.name

sc_folder = glob.glob(folder + "/run_*/busco_sequences/single_copy_busco_sequences")[0]
mc_folder = glob.glob(folder + "/run_*/busco_sequences/multi_copy_busco_sequences")[0]
frag_folder = glob.glob(folder + "/run_*/busco_sequences/fragmented_busco_sequences")[0]

tmp_count = glob.glob(sc_folder + "/*.fna")

if len(tmp_count)==0:
	if args.assembly==None:
		print("No assembly specified, but is required because no fna files exist")
		print("Use -a/--assembly flag to specify input assembly")
	else:
		assembly_name = os.path.abspath(args.assembly)

delim = args.delim
cpus = args.cpu

print("""

Clean BUSCO results
1. Check multi-copy busco sequences for 100% duplicates
	- This will result in more single-copy genes
2. Rename sequences and concatenate

Results will have the format:
>locus_name|sample_name

Delimiter | can be changed with the -d flag.

""")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting BuscoCleaner...")
print("\tBUSCO output folder -> "+ args.folder)
print("\tSample name -> "+ args.name)
print("\tOriginal assembly -> "+ args.assembly)
print("\tDelimiter -> " + delim)
print("\tNumber of CPU -> "+ str(cpus))

########################################
################# CODE #################
########################################

#### CREATING FNA FILES
if len(tmp_count)==0:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Creating fna files :::")
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Reading Original Assembly :::")
	
	assembly=list(SeqIO.parse(assembly_name,"fasta"))

	sc_buscos = glob.glob(sc_folder + "/*.faa")
	mc_buscos = glob.glob(mc_folder + "/*.faa")
	frag_buscos = glob.glob(frag_folder + "/*.faa")
	
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Grabbing Single-Copy Loci :::")
	p_map(Distributor, sc_buscos, num_cpus=cpus)
	gc.collect()
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Grabbing Multi-Copy Loci :::")
	p_map(Distributor, mc_buscos, num_cpus=cpus)
	gc.collect()
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Grabbing Fragmented Loci :::")
	p_map(Distributor, frag_buscos, num_cpus=cpus)
	del assembly
	gc.collect()

#### CHECKING FOR DUPLICATES IN MULTI-COPY FNA FILES
print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Checking for 100% Duplicates in Multi-Copy Sequence Folder :::")
mc_fnas1 = glob.glob(mc_folder+"/*.fna")
sc_fnas1 = glob.glob(sc_folder+"/*.fna")

mkdir_p(mc_folder + "/Duplicates")
results=p_map(mcChecker, mc_fnas1, num_cpus=cpus)

sc_fnas2 = glob.glob(sc_folder+"/*.fna")
mc_fnas2 = glob.glob(mc_folder+"/*.fna")
dup_count = glob.glob(mc_folder+"/Duplicates/*.fna")

print("Single Copy Busco Loci: "+ str(len(sc_fnas1)))
print("Multi Copy Busco Loci: "+ str(len(mc_fnas1)))
print("::::::::::::::")
print("Duplicates in Multi Copy Busco Loci: "+ str(len(dup_count)))
print("::::::::::::::")
print("New Single Copy Busco Loci: "+ str(len(sc_fnas2)))
print("New Multi Copy Busco Loci: "+ str(len(mc_fnas2)))


#### RENAMING AND CONCATENATING
print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Renaming and Concatenating Loci into " + folder + "/BuscoCleaner :::")
mkdir_p(folder+"/BuscoCleaner")
results=p_map(Renamer, sc_fnas2, num_cpus=cpus)

if args.output!=None:
	sp.call("cat " + folder + "/" + name + "+_busco.fasta >> " + os.path.join(os.path.abspath(args.output),"busco_loci.fasta"), shell=True)