#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 0-01:00:00
#SBATCH --array 1-21
#SBATCH -J msprime_sim
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# wrapper for the Ychromosome simulations
outfile=$1
if (( ${SLURM_ARRAY_TASK_ID} == 1 ))
then
	echo -e "Count\tReplicate\tGeneFlowProp" > ${outfile}
fi

# fetch the gene flow proportion to test from file
i=$(cat array_proportions.txt | grep -v "#" | awk -v line=${SLURM_ARRAY_TASK_ID} 'NR == line {print $2}')

for REP in $(seq 1 100)
do
python3 denti_simulations_ychrom.py ${i} out/ychrom_sims_rep_${REP}_prop_${i}
# parse the output to a master file
## fetch the counts of monophyletic denti/mitis trees, if any
COUNT=$(cat out/ychrom_sims_rep_${REP}_prop_${i}_distance_counts.txt | awk ' $1 == "2" ' | cut -f 2)
# if count variable is empty, set it to 0
if [ -z "${COUNT}" ]
then
	COUNT=0
fi
# print the distance, count, replicate, and geneflow proportion to the master file
echo -e "${COUNT}\t${REP}\t${i}" >> ${outfile}
done

