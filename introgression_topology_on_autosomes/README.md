This directory contains the scripts used for identifying the Y-tree topology (denti/mitis monophyly) in the autosomes.

## Dependencies

The scripts uses the following software/packages (versions used in parenthesis):
- PhyML (3.3.20190321)
- the python module ete3 (3.1.2)
- bcftools (1.20)
- My phylogenomics repo: https://github.com/axeljen/phylogenomics. After cloning and installing dependencies, the path to this repo should be modified in the codeml.py script.

## Workflow

1. Build phylogenies in 10 kb, nonoverlapping, sliding windows along the genome.

```bash

# path to directory with chromosome specific vcf files
vcfdir=../vcffiles/
# file with samples to include
samples=samples.txt
# output directory
outdir=10K_nj_trees
mkdir -p ${outdir}

# run it for each chrom
for vcf in $(find ${vcfdir} -name "*.vcf.gz")
do
chrom=$(bcftools view -H ${vcf} | head -n 1 | cut -f 1)
sbatch -J ${chrom}.phyml phyml_nj_windows.sh \
	${vcf} \
	${outdir}/${chrom}.txt \
	${samples}
done

```

2. Identify trees where denti and mitis are monophyletic, nested within a monophyletic mitis clade, to the exclusion of a monophyletic mona clade. (the other taxa are not considered in this analysis)

```bash
# loop through trees from all chromosomes
for tree in $(find ${outdir} -name "*.txt" | grep -v "log")
do
echo ${tree}
python3 find_mitis_to_denti_intro_trees.py --treefile ${tree} --samples popfile_trees.txt -o denti_mitis_monophyletic_trees.txt
done
```
