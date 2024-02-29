#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH -J bpp
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# give controlfile as command line input
controlfile=${1}

cd $(dirname ${controlfile})

bpp --cfile $(basename ${1})