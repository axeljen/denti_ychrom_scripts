#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH -J iqtree
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

set -e

#Specify how many threads we're using
THREADS=16

# load modules
module load gcc openmpi
module load bioinfo-tools iqtree python/3.8.7

# take a set of parameters as command line inputs
VCF=$1
OUTPUT=$2
SAMPLES=$3

# make dir of output if it doesn't exist
mkdir -p $(dirname ${OUTPUT})

# change bootstraps and model here if wanter
BOOTSTRAPS=1000
MODEL="GTR"
WINDOW=25000
STEP=500000
MAX_MISSING=0.1

# path to python script
PYTHONSCRIPTS=/home/axeljen/phylogenomics

# start running
python3 ${PYTHONSCRIPTS}/windowPhylogenies.py -i ${VCF} \
	-o ${OUTPUT} \
	--sequences ${SAMPLES} \
	--max-missing ${MAX_MISSING} \
	--bootstraps ${BOOTSTRAPS} \
	--model ${MODEL} \
	-w ${WINDOW} \
	--step-size ${STEP} \
	-T 20 \
	--haploidize \
	--heterozygotes IUPAC