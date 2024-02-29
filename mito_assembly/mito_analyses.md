# Mitochodnrial assembly and phylogeny

This directory contains the scripts used for assemblying, annotating and analyzing the mitochondrial genomes.

I had the reads stored in .cram format (mapped to reference), so that's the starting point here. From there, the sripts were used in the following order:

1. mitofinder_asm.sh
	This makes the first, "raw" assembly of the mtDNA for each sample.
2. rotate_and_align_mt_genomes.sh
	This script contains some steps for processing and making a first alignment of the mtDNA genomes. This alignment was then used to manually check the integrity of each genome, and I subsequently extracted all individual sequences again for annotation.
3. mitofinder_annotate.sh
	This script annotates the rotated and trimmed mtDNA genomes.
4. extract_mt_features.sh 
	Contains commands for extracting the annotated features from the mtDNA genomes.
5. align_and_partition_features.sh
	Aligns these features first individually, and then concatenates the alignments into a single file, accompanied by a partitionmap file for use in phylogenetic analysis.
6. iqtree_partitioned_phylogeny.sh
	Runs the phylogenetic analysis in iqtree, with the builtin function modelfinder to find the best model for each partition.