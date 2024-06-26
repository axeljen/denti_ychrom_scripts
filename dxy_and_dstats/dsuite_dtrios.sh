#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-05:00:00
#SBATCH -J dsuite
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

module load python3
module load bioinfo-tools bcftools

INPUT_VCF=$1

SETS=$2

OUT=dstat-all-samples

mkdir -p ${OUT}

TREE=all_samples_topology.nwk

CHROM=$(bcftools view -H ${INPUT_VCF} | head -n 1 | cut -f 1)

Dsuite Dtrios -o ${OUT}/$CHROM -t ${TREE} $INPUT_VCF $SETS

