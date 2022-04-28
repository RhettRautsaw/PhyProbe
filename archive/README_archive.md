

1. Assembly: rnaSpades or Trinity
2. BUSCO 4.1.4: tetrapoda_odb10
	- BuscoCleaner.py:
		- Check for duplicates in multi-copy and move to single copy
		- Concatenate all individuals into busco_loci.fasta
3. CodAn: VERT_full
4. Extract longest isoforms (LongForms.py)
5. cd-hit-est: 100% 
6. Extract sqCL loci (LocusExtractor.py)
7. Orthofinder on clustered ORFs
	- scOrthoSorter.py:
		- Extract single copy loci with less than 50% missing individuals
8. Combine Data Sources
	- UniqRenamer.py:
		- Compare BUSCO loci to sqCL targets.
		- Extract matching regions and rename appropriately
		- This makes sure that we don't have repeat loci for phylogenetics later
	- UniqRenamer.py:
		- Compare Orthofinder loci to combined sqCL and BUSCO loci
	- SampleSorter.py:
		- Sort combined BUSCO, OF, and sqCL extracted loci into each individual
9. Phasing (Phaser.py)
