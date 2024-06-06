This directory contains the scripts and files used to calculate nucleotide divergence (dxy) and D-statistics.

## Dependencies
These scripts use the following software/external scripts (versions used in parentheses):
- pixy (1.2.5.beta1)
- bcftools (1.20)
- Dsuite (0.4r42)
- ABBABABAwindows.py from Simon Martin's genomics_general repo (https://github.com/simonhmartin/genomics_general)

Dxy was calculated between male samples among the mona and mitis groups, using the pixy_dxy.sh script together with the pops-denti-mitis-dxy-males.txt popfile:


```bash

# path to directory with chromosome VCFs
vcfdir=../vcffiles
# popfile
popfile=pops-denti-mitis-dxy-males.txt
# window size
window=10000
# output name
outname=dxy
# run pixy
for vcf in $vcfdir/*.vcf.gz
do
sbatch pixy_dxy.sh $vcf $popfile $window $outname
done

```

D-stats was first calculated on a chrom-by-chrom basis with dsuite_dtrios.sh, and the chromosome-specific outputs were combined with dtrios_combine.txt. The file sets_samples.txt was used as the setsfile, and the all_samples_topology.nwk as the treefile.

```bash

# path to vcf files (note that this should be run on biallelic SNPs only)
vcfdir=../vcffiles
# sets file
setsfile=sets_samples.txt
# topology file
treefile=$(realpath all_samples_topology.nwk)

# first run Dsuite per chromosome
for vcf in $vcfdir/*.vcf.gz
do
sbatch dsuite_dtrios.sh $vcf $setsfile
done

# when these are finished, use dtrios combine to combine the dstatistcs from all chromosomes
files=$(find dstat-all-samples -name "*Dmin.txt" | rev | cut -d "_" -f2- | rev)
Dsuite DtriosCombine --out-prefix dstat-all-samples/dstat-allsamples-combined -t ${treefile} ${files}

```


The fd_windows.sh script was used to calculate Fd/dstatistics in 10k sliding windows along the genome, together with the fd_popfile.txt.

```bash
# path to vcf files (note that this should be run on biallelic SNPs only)
vcfdir=../vcffiles
# popfile and population/trio set up, window size etc is specified in the fd_windows.sh script,
# so just run the script for each vcf
for vcf in $vcfdir/*.vcf.gz
do
sbatch fd_windows.sh $vcf
done

```
