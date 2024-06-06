#!/bin/bash

#SBATCH -A naiss2023-5-506
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
module load bioinfo-tools PhyML python/3.8.7

# take a set of parameters as command line inputs
VCF=$1
OUTPUT=$2
SAMPLES=$3
#OUTGROUP=$4
#WINDOWS=$5
# make dir of output if it doesn't exist
mkdir -p $(dirname ${OUTPUT})

# change bootstraps and model here if wanter
MODEL="GTR"
WINDOW=10000
STEP=10000
MAX_MISSING=0.1
BOOTSTRAPS=0
# path to python script
PYTHONSCRIPTS=/home/axeljen/phylogenomics

# start running
python3 ${PYTHONSCRIPTS}/windowPhylogenies.py -i ${VCF} \
	-o ${OUTPUT} \
	--phyml \
	--optimize n \
	--sequences ${SAMPLES} \
	--max-missing ${MAX_MISSING} \
	--bootstraps ${BOOTSTRAPS} \
	--model ${MODEL} \
	-w ${WINDOW} \
	--step-size ${STEP} \
	-T 20 \
	--haploidize \
	--heterozygotes IUPAC