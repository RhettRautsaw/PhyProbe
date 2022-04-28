# PhyProbe Scripts

This folder has all the scripts used by `PhyProbe`. However, several of the other scripts that PhyProbe calls are useful on their own. 

`PhyProbe.py`: The main script used for the full pipeline on the previous page. 

<br>

***

| Script				| Description                                                              |
|-----------------------|--------------------------------------------------------------------------|
|`AliBaSeqCleaner.py`	| Cleans the results of AliBaSeq in several ways: <br/> (1) Remove sequences with coverage less than a certain % of the reference locus. <br/> (2) Add sample name to fasta headers. <br/> (3) Concatenate loci into one sample fasta |
|`BuscoCleaner.py`		| Cleans the results of BUSCO in several ways: <br/> (1) If the fna files were not produced, it will use the information from the faa files and the original assembly to extract them. <br/> (2) Check the multi-copy sequence folder for sequences with 100% identical duplicates (not true multi-copy loci...only assembly artifacts) <br/> (4) Add sample name to fasta headers. <br/> (5) Concatenate loci into one sample fasta |
| `DropGappy.py`		| Calculate % of alignment that is made up of gaps for each sample and remove if greater than a given threshold |
| `FastaRenamer.py`		| Rename fasta headers numerically with given prefix (e.g., contig0, contig1, contig2...) |
| `LongForms.py`		| Extract the longest isoform from Trinity or rnaSpades assembly |
| `MCL.py`				| Perform Markov Clustering on a provided fasta via self-BLAST search |
| `MissingDataFilter.py`| Subset alignments for a specified amount of missingness on two different axes: <br/> (1) Filter for genes with a certain % of total taxa <br/> (2) Filter for taxa with a certain % of total genes |
| `OrthoCleaner.py`		| Cleans the results of Orthofinder to extract single-copy orthologous sequences. Orthofinder does this, but requires 100% taxa coverage. OrthoCleaner will allow identification of loci that are present in a specified % of taxa |
| `SeqSorter.py`		| Separate a fasta into multiple fastas based on specified delimiter in the fasta header and position. For example, with fasta header `>gene\|sample`, SeqSorter can create separate fastas for each gene present with each containing all the samples (useful for phylogenetics). Alternatively, you can create separate fastas for each sample with each containing all the genes for that sample. |
| `Ungap.py`			| Remove gaps from a fasta (unalign) |
