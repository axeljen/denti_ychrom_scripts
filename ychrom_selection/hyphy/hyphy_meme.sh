#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 0-10:00:00
#SBATCH -J hyphy_meme
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# load modules
module load bioinfo-tools HyPhy

# alignment to test
ALN=ygenes_concatenated.phy

# annotated species tree
TREE=hyphy_tree.nwk

hyphy meme --alignment ${ALN} --tree ${TREE} CPU=20