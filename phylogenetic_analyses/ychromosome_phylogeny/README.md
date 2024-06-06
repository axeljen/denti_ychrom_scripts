# Y chromosome phylogeny

The file ychrom_aln.phy contans a multiple sequence alignment of the Y chromosome including all males, contructed as consensus sequences from the genotyped called against the rhesus macaque reference. Sites with more than 10 % missing data were excluded.

Build a phylogeny with the iqtree.sh script:

```bash
# make a directory to store output files
mkdir -p out
# run script
sbatch iqtree.sh ychrom_aln.phy out/ychrom_phylogeny

```