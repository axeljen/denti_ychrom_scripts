#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2-00:00:00
#SBATCH -J codeml
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

set -e

# load modules
module load bioinfo-tools paml

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