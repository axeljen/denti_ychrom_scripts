#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J slim
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# load modules
module load bioinfo-tools SLiM/4.0.1


# script with model
SLIMSCRIPT=${1}

# variables from command line
SEL=${2}
POPSIZE=${3}
INTROPROP=${4}
OUTDIR=${5}

mkdir -p ${OUTDIR}



# running reps sequenctailly
NREPS=2


for i in $(seq 1 ${NREPS})
do
name=${SLURM_JOB_ID}

# run slim
slim \
	-d name=${name} \
	-d iteration=${i} \
	-d s=${SEL} \
	-d popsize=${POPSIZE} \
	-d introprop=${INTROPROP} \
	${SLIMSCRIPT}

# move output to results folder
mv ${name}_${i}_N-${POPSIZE}_Intro-${INTROPROP}_s-${SEL}.txt ${OUTDIR}

done
