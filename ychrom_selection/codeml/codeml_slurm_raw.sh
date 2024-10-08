#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 0-05:00:00
#SBATCH -J codeml
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

set -e

# load modules
module load bioinfo-tools paml/4.9j

# take some arguments from the command line
aln=$1
tree=$2
samples=$3
name=$4
outgroup=$5

foreground="FOREGROUND_BRANCHES"
echo ${foreground}

# make the output directory
outdir=out/${name}

mkdir -p ${outdir}

python3 codeml.py \
	--alignment $aln \
	--tree $tree \
	--samples $samples \
	-w $outdir \
	--name $name \
	-f ${foreground} \
	--outgroup ${outgroup} \
	--maxmissing 0.5