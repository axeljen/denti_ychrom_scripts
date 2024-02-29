#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 0-12:00:00
#SBATCH -J astral
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/astral-%x-%j.out
#SBATCH -e ./logs/astral-%x-%j.error

module load bioinfo-tools java

ASTRAL_PATH=/home/axeljen/ASTRAL/ASTRAL/Astral
# input file with trees as first cmnd input
INPUT=$1
#Give output directory to store trees in as second cmd input
OUTPUT=$2
#change node annotations here
ANNOTATION=$3

java -jar $ASTRAL_PATH/astral.5.7.4.jar -t $ANNOTATION -i $INPUT -o ${OUTPUT}