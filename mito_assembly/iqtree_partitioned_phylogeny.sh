#!/bin/sh

#SBATCH -A snic2022-5-561
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

# samples to include in the analysis
SAMPLES=$1

# name of run
RUN_NAME=$2

# path to alignment
RAW_ALN=/crex/proj/sllstore2017021/nobackup/DENTI_LESULA_PROJECT/ANALYSES/MITOGENOMES_NEW/TRIMMED_POLISHED_ALIGNMENTS/trimmed_concatenated_features_aligned.phy

# path to partitions
PARTITIONS=/crex/proj/sllstore2017021/nobackup/DENTI_LESULA_PROJECT/ANALYSES/MITOGENOMES_NEW/TRIMMED_POLISHED_ALIGNMENTS/trimmed_concatenated_features_aligned_raxml_partitions.txt

# make directory for run
mkdir -p ${RUN_NAME}

# first subset alignment
python3 ~/phylogenomics/subsetSequence.py -i ${RAW_ALN} -s ${SAMPLES} -o ${RUN_NAME}/${RUN_NAME}_sub_align.phy

# then do a modeltest
iqtree2 -s ${RUN_NAME}/${RUN_NAME}_sub_align.phy -p ${PARTITIONS} -m TESTMERGEONLY --prefix ${RUN_NAME}/${RUN_NAME}_modeltest

# then do a tree search using the best scheme from modeltest
iqtree2 -s ${RUN_NAME}/${RUN_NAME}_sub_align.phy -p ${RUN_NAME}/${RUN_NAME}_modeltest.best_scheme.nex --prefix ${RUN_NAME}/${RUN_NAME}_tree -B 1000