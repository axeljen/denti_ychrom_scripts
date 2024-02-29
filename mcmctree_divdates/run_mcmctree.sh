#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5-00:00:00
#SBATCH -J mcmctree
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/mcmctree-%x-%j.out
#SBATCH -e ./logs/mcmctree-%x-%j.error

module load bioinfo-tools paml/4.9j

#give controlfile as input
CONFIG_FILE=$1

cd $(dirname ${CONFIG_FILE})

mcmctree $(basename ${CONFIG_FILE})
