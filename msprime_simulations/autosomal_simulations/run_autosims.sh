#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 0-01:00:00
#SBATCH -J jobname
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

module load bioinfo-tools

# this script will simulate 100 "chromosomes" of 1 Mb each 
# then, we'll calculate Dstatistics for mitis -> denti on each chromosome separately, and then combine them
# this will yield one replicate simulation of the Dstatistics for each gene flow proportion

# specify the gene flow proportion to test
i=$1

# output name as second, this is the filename in to which the output will be concatenated
OUTNAME=$2

# number of "chromosomes" to simulate
NREPS=100

# make a directory for the vcf files of this repetition
mkdir -p vcf_files/${SLURM_JOB_ID}

# simulate 100 chromosomes
python3 denti_simulations_autosomes.py ${i} vcf_files/${SLURM_JOB_ID}/sims_rep_${SLURM_JOB_ID}_prop_${i} ${NREPS}

### Dstats section

# make a directory for the dstat files of this repetition
mkdir -p dstatfiles/${SLURM_JOB_ID}

for REP in $(seq 1 ${NREPS})
do
Dsuite Dtrios -o dstatfiles/${SLURM_JOB_ID}/dstats_${SLURM_JOB_ID}_${REP}_prop_${i} -t Dstat_tree.txt vcf_files/${SLURM_JOB_ID}/sims_rep_${SLURM_JOB_ID}_prop_${i}_rep${REP}.vcf Dstat_SETS.txt
done

# find all the combine files 
for REP in $(seq 1 ${NREPS})
do
find dstatfiles/${SLURM_JOB_ID} -name "*_Dmin.txt" | rev | cut -d "_" -f 2- | rev > dstatfiles/${SLURM_JOB_ID}/Dmin_files_${SLURM_JOB_ID}_prop_${i}.txt
done

DMINFILES=$(cat dstatfiles/${SLURM_JOB_ID}/Dmin_files_${SLURM_JOB_ID}_prop_${i}.txt)

Dsuite DtriosCombine -t Dstat_tree.txt -o dstatfiles/${SLURM_JOB_ID}/combined_Dmin_${SLURM_JOB_ID}_prop_${i} ${DMINFILES}

# parse the output to a master file
python3 parse_dsuite_output.py dstatfiles/${SLURM_JOB_ID}/combined_Dmin_${SLURM_JOB_ID}_prop_${i}_combined_tree.txt ${SLURM_JOB_ID} ${i} >> ${OUTNAME}

# remove the vcf files
rm -r vcf_files/${SLURM_JOB_ID}