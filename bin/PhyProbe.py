#!/usr/bin/env python3

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import multiprocessing as mp
from psutil import virtual_memory
import glob

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

       __   ____  _           ____            _           
   ___|    |  _ \| |__  _   _|  _ \ _ __ ___ | |__   ___  
  |   |__  | |_) | '_ \| | | | |_) | '__/ _ \| '_ \ / _ \ 
--|        |  __/| | | | |_| |  __/| | | (_) | |_) |  __/ 
  |______  |_|   |_| |_|\__, |_|   |_|  \___/|_.__/ \___| 
                        |___/                             

RNA-Seq (Individual Samples):
	1. Transcriptome Assembly (rnaSpades)
	2. BUSCO & BuscoCleaner.py
	3. CDS Prediction (CodAn)
	4. Extract Longest Isoforms (Isoformer.py)
	5. Cluster Sequences (cd-hit-est)
	6. Extract Target Loci (LocusExtractor.py)

RNA-Seq (Grouped Samples):
	NOTE: This requires that steps 3-5 are complete for all RNA-Seq
	WARNING: This will likely be computationally intensive if you have a lot of samples
	7. Orthofinder & scOrthoSorter.py
	8. Combine Data Sources (BUSCO, Orthofinder, and LocusExtractor)

RNA-Seq (Individual Samples):
	9. Phasing (Phaser.py)

Seq-Cap and Other Datasources (Individual Samples):
	Combine original targets with newly idenfied RNA-Seq loci, then do steps:
	6. Extract Target Loci (LocusExtractor.py)
	9. Phasing (Phaser.py)

Summarize:
	10. Combine Data Sources (Combiner2.py)
	11. Align and Infer Phylogeny (Phylogenify.py)
		- MAFFT
		- Trimal
		- SpruceUp
		- IQTREE
		- TreeShrinks
	12. ASTRAL-MP Species Tree Inference

""")

parser.add_argument("-i","--input",
					type=str,
					help="Folder with DNA/RNA read data and/or previously completed steps. Folder name should correspond to sample name and will be used as a prefix. \n If running step 7, this should be a folder containing the clustered fastas from all RNA-Seq samples.")
parser.add_argument("-s", "--steps",
					type=str,
					help="Comma-separated list of steps in the pipeline to perform. (default: %(default)s) \n Note that each step requires that files and folder structure be correct.")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Delimiter used to separate gene or locus names from sample names (default: %(default)s)")
parser.add_argument("-c","--cpu",
					type=int,
					default=mp.cpu_count(),
					help="Number of threads to be used in each step. (default: All cpus = %(default)s)")
parser.add_argument("-mem","--memory",
					type=int,
					default=(round(virtual_memory().total/1e9)-1),
					help="Required for Step 1. Max amount of memory for assembly. (default: All mem = %(default)s)")
parser.add_argument("-l","--lineage",
					type=str,
					help="Required for Step 2. Path to BUSCO lineage or name of BUSCO lineage")
parser.add_argument("-m","--model",
					type=str,
					help="Required for Step 3. Path to CodAn model for transcriptome CDS prediction")
parser.add_argument("-t","--targets",
					type=str,
					help="Required for Steps 7 & 8. Fasta with targets from which to extract loci")
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

if args.input==None:
	print("No input specified. Use -i/--input flag to specify input")
else:
	input = os.path.abspath(args.input)
	prefix = args.input

if args.steps==None:
	print("No steps specified. Use -s/--steps flag to specify steps")
	quit()
else:
	steps = args.steps.split(",")

delim = args.delim
cpus = args.cpu
mem = args.memory

start_dir=os.getcwd()

print("""

       __   ____  _           ____            _           
   ___|    |  _ \| |__  _   _|  _ \ _ __ ___ | |__   ___  
  |   |__  | |_) | '_ \| | | | |_) | '__/ _ \| '_ \ / _ \ 
--|        |  __/| | | | |_| |  __/| | | (_) | |_) |  __/ 
  |______  |_|   |_| |_|\__, |_|   |_|  \___/|_.__/ \___| 
                        |___/                             

""")

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting PhyProbe v1.0 :::")
print("\tInput -> "+ args.input)
print("\tSteps -> "+ args.steps)
print("\tDelimiter -> " + delim)
print("\tNumber of CPU -> "+ str(cpus))

if '1' in steps:
	print("\tMax Memory -> "+ str(mem))
if '2' in steps:
	if args.lineage==None:
		print("No BUSCO lineage specified. Using auto-mode...this might take a while.")
		print("To avoid, download desired database and use -l/--lineage flag to specify input")
		lineage = 'auto'
	else:
		lineage = os.path.abspath(args.lineage)
	print("\tLineage -> "+ lineage)
if '3' in steps:
	if args.model==None:
		print("No CodAn model specified. Use -m/--model flag to specify input")
		quit()
	else:
		model = os.path.abspath(args.model)
		print("\tCodAn Model -> "+ model)
if '6' in steps or '8' in steps:
	if args.targets==None:
		print("No target fasta specified. Use -t/--targets flag to specify input")
		quit()
	else:
		targets_name = os.path.abspath(args.targets)
		print("\tTargets Fasta -> "+ args.targets)

########################################
################# CODE #################
########################################

os.chdir(input)
if '1' in steps or '6' in steps or '9' in steps:
	reads = [os.path.abspath(x) for x in glob.glob("*.fastq*")+ glob.glob("*.fq*")]
	reads.sort()

### STEP 1 : ASSEMBLY
if '1' in steps:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Step 1: rnaSpades Assembly :::')
	if len(reads)==1:
		sp.call("rnaspades.py -s " + reads[0] + " -t " + str(cpus) + " -m " + str(mem) + " -o 01_assembly", shell=True)
	elif len(reads)==2:
		sp.call("rnaspades.py -1 " + reads[0] + " -2 " + reads[1] + " -t " + str(cpus) + " -m " + str(mem) + " -o 01_assembly", shell=True)
	elif len(reads)>2 or len(reads)==0:
		print(str(len(reads)) + " fastq files identified")
		print("Please only include one file for single-end mode or two files for paired-end mode")
		quit()

### STEP 2 : BUSCO
if '2' in steps:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Step 2: BUSCO Prediction :::')
	sp.call("busco -i 01_assembly/transcripts.fasta -c " + str(cpus) + " -o 02_busco -m geno -l " + lineage, shell=True)
	sp.call("BuscoCleaner.py -f 02_busco -n " + prefix + " -c " + str(cpus), shell=True)
	sp.call("cat 02_busco/" + prefix + ".buscos.fasta >> " + start_dir + "/busco_loci.fasta", shell=True)

### STEP 3 : CODAN
if '3' in steps:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Step 3: CodAn Prediction :::')
	sp.call("codan.py -t 01_assembly/transcripts.fasta -m " + model + " -o  03_codan -c " + str(cpus), shell=True)

### STEP 4 : ISOFORMER
if '4' in steps:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Step 4: Extracting Longest Isoforms :::')
	mkdir_p("04_long_isoforms")
	sp.call("LongForms.py -i 03_codan/ORF_sequences.fasta -o 04_long_isoforms --auto", shell=True)

### STEP 5 : CD-HIT-EST
if '5' in steps:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Step 5: Clustering Isoforms at 100% Identity :::')
	mkdir_p("05_clustered")
	sp.call("cd-hit-est -i 04_long_isoforms/ORF_sequences.long.fasta -o 05_clustered/" + prefix + ".clust.fasta -c 1 -d 0", shell=True)
	sp.call("perl -pi -e 's/>/>" + prefix + delim +"/g' 05_clustered/" + prefix + ".clust.fasta", shell=True)

### STEP 6 : LOCUS-EXTRACTOR
if '6' in steps:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Step 6: Locus Extractor :::')
	sp.call("LocusExtractor.py -r " + ' '.join(reads) + " -t " + targets_name + " -d \"" + delim + "\" -p " + prefix + " -o 06_targets -c " + str(cpus), shell=True)
	sp.call("cat 06_targets/*.targets.fasta >> " + start_dir + "/target_loci.fasta", shell=True)

### STEP 7 : ORTHOFINDER
if '7' in steps:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Step 7: Orthofinder :::')
	os.chdir(start_dir)
	sp.call('orthofinder -f ' + input + '-d -og -o 06_orthofinder -t ' + str(cpus), shell=True)
	sp.call('scOrthoSorter.py -f 06_orthofinder -c ' + str(cpus), shell=True)
	sp.call('cp 07_orthofinder/scOS_combined.fasta ' + start_dir + '/orthofinder_loci.fasta', shell=True)

### STEP 8: DATA COMBINER
if '8' in steps:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Step 8: Combine data sources :::')
	os.chdir(start_dir)
	sp.call('UniqRenamer.py -i busco_loci.fasta -t ' + targets_name + ' -c ' + str(cpus), shell=True)
	sp.call('cat ' + targets_name + ' busco_loci_renamed.fasta > combined_loci.fasta', shell=True)
	sp.call('UniqRenamer.py -i orthofinder_loci.fasta -t combined_loci.fasta -c ' + str(cpus), shell=True)
	sp.call('cat combined_loci.fasta orthofinder_loci_renamed.fasta target_loci.fasta > combined_loci2.fasta', shell=True)
	sp.call('SeqSorter.py -i combined_loci2.fasta -c ' + str(cpus), shell=True)
	os.chdir("SeqSorter")
	files = glob.glob("*.fasta")
	def LongForms(file):
		print(file)
		sp.call("LongForms.py -i " + file, shell=True)
	
	results = p_map(LongForms, files, num_cpus=cpus)


### STEP 9: PHASING
if '9' in steps:
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Step 9: Phase Loci :::')
	sp.call("Phaser.py -r " + ' '.join(reads) + " -t " + targets_name + " -o 09_phasing -c " + str(cpus), shell=True)



#mafft --auto --adjustdirectionaccurately --thread 16 {} > {}.aln
#trimal -in {}.aln -out {}.trim -automated1
#sed 's/PLACEHOLDER/{}.trim/g' config_example.conf > {}.spruceup.conf
#python -m spruceup {}.spruceup.conf
#mv 0.99_lognorms-cutoff-{}.spruceup.fasta {}.spruceup.fasta