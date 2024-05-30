
# Divergence date estimates with mcmctree

This directory contains scripts/controlfiles used for the mcmctree analyses of autosomal and Y-chromosomal data.

Alignments were constructed from the vcf-files in the same way as for the bpp analyses (but only 10 loci á 10 kb were sampled), and then two mcmctree-runs per analysis were run, using the run_mcmctree.sh and *_controlfile.ctl files. For the autosomes, the loci sampling and mcmctree analysis was repeated 10 times, and the runs were merged with mcmctree by changing the "print" parameter to "-1".

The script is executed as:
	sbatch run_mcmctree.sh controlfile.ctl