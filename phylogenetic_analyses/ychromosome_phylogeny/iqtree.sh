#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 1-00:00:00
#SBATCH -J iqtree
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# load modules
module load gcc openmpi
module load bioinfo-tools iqtree python/3.8.7

#input fasta alignment
INPUT=$1

#output prefix
OUTPUT=$2

#nThreads
THREADS=20

iqtree2 -s $INPUT -m MFP -pre $OUTPUT -nt $THREADS -B 1000
