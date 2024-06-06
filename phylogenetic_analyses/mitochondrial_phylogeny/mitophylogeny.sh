#!/bin/sh

#SBATCH -A naiss2023-5-506
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 0-10:00:00
#SBATCH -J iqtree
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# using iqtree to find the best partitioning scheme for all my partitions
module load bioinfo-tools iqtree

# name of run
RUN_NAME=$1

# path to alignment
RAW_ALN=alignment/trimmed_concatenated_features_aligned.phy

# path to partitions
PARTITIONS=alignment/trimmed_concatenated_features_aligned_raxml_partitions.txt

# make directory for run
mkdir -p ${RUN_NAME}

# run modeltest
iqtree2 -s ${RAW_ALN} -p ${PARTITIONS} -m TESTMERGEONLY --prefix ${RUN_NAME}/${RUN_NAME}_modeltest

# then do a tree search using the best scheme from modeltest
iqtree2 -s ${RAW_ALN} -p ${RUN_NAME}/${RUN_NAME}_modeltest.best_scheme.nex --prefix ${RUN_NAME}/${RUN_NAME}_tree -B 1000