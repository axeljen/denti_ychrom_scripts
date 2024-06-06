# Autosomal phylogeny scripts

This directory contains the scripts to build the astral phylogeny (Figure 1A).

## Dependencies

- astral (https://github.com/smirarab/ASTRAL, I used v. 5.7.4.)
- IQtree (http://www.iqtree.org/, I used v. 2.2.2.6-omp-mpi)

- My phylogenomics repo (https://github.com/axeljen/phylogenomics). After cloning this, adjust the PYTHONSCRIPTS path in iqtree_sliding_windows.sh to point to the correct location.

## Workflow

1. Run iqtree_sliding_windows.sh to build the gene/window trees.
I ran it on chromosome-wise vcf files with genotypes (SNPs and invariants):

```bash
	VCF_DIR=../../vcffiles
	OUT_DIR=out
	samples=samples.txt # this is just a file listing all the samples included in the vcf file (or a subset if wanted)
	for VCF in $VCF_DIR/*.vcf.gz; do
		# using bcftools to get the chromosome name from the vcf file, requires bcftools in PATH 
		chrom=$(bcftools view -H ${VCF} | head -n 1 | cut -f 1)
		sbatch iqtree_sliding_windows.sh ${VCF} ${OUT_DIR}/${chrom}_trees.txt ${samples}
	done
```
2. when all jobs are done, merge the trees across all chromosomes:
```bash
	head -n 1 $(find ${OUT_DIR} -name "*trees.txt" | head -n 1) > ${OUT_DIR}/all_trees.txt
	for f in $(find ${OUT_DIR} -name "*trees.txt" | grep -v all_trees.txt); do
		# remove NA tree and concatenate the rest to the joint file
		tail -n +2 $f | awk ' $7 != "NA" ' >> ${OUT_DIR}/all_trees.txt
	done
```
3. Now run astral on the merged trees:
```bash
	# get only the tree column without header to a new file
	cat ${OUT_DIR}/all_trees.txt | tail -n +2 | cut -f 7 > ${OUT_DIR}/all_trees.trees
	# run astral with LPP support annotations
	sbatch run_astral.sh ${OUT_DIR}/all_trees.trees ${OUT_DIR}/astral_output.tre 3

```