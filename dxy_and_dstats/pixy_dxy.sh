#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 1-00:00:00
#SBATCH -J pixy_dxy
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# load modules
module load bioinfo-tools pixy bcftools

# how many threads are we going at?
THREADS=20

# vcf file as first input
VCF=$1

# popfile as second
POPS=$2

# window size as third
WSIZE=$3

# filename as fourth
NAME=$4

# grab chromname from first vcf record
chrom=$(bcftools view -H ${VCF} | head -n 1 | cut -f 1)

# make an output dir from this
mkdir -p ${NAME}_${WSIZE}-windows

# and run pixy
pixy --stats dxy \
	--vcf ${VCF} \
	--populations ${POPS} \
	--output_folder ${NAME}_${WSIZE}-windows \
	--output_prefix ${chrom} \
	--n_cores ${THREADS} \
	--window_size ${WSIZE}

