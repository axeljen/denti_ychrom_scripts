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
NREPS=100


# make a unique name of simulation from slurm job id
name=${SLURM_JOB_ID}

for i in $(seq 1 ${NREPS})
do

# run slim
slim \
	-d name="${name}" \
	-d iteration=${i} \
	-d s=${SEL} \
	-d popsize=${POPSIZE} \
	-d introprop=${INTROPROP} \
	${SLIMSCRIPT}

done

# move all output to results folder
mv ${name}_*.txt ${OUTDIR}
