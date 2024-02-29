# Divergence date estimates with mcmctree

This directory contains scripts/controlfiles used for the mcmctree analyses of autosomal and Y-chromosomal data.

Alignments were constructed from the vcf-files, and then two mcmctree-runs per analysis were run, using the run_mcmctree.sh and *_controlfile.ctl files.

For the autosomal data, I ran this script for ten different data sets containing ten loci each. These runs were then merged using "print = -1".