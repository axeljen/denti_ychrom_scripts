#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 0-00:30:00
#SBATCH -J mitofinder
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

#mitofinder is installed manually, specify the path to your installation here NOTE: you need to change this to your own installation
MITOFINDER_PATH=/domus/h1/axeljen/MitoFinder

#specify the path to the reference genome to work on, instructions on how to find and download this is available on mitofinder's github
REFERENCE=/crex/proj/sllstore2017021/nobackup/GUENON_GENOMES/REFERENCE/MITOGENOME/MITOFINDER_REFERENCE/chl.sab_NC_008066.1.mito.gb

#give mitochondrial fasta sequence as input
MITOGENOME=$1

FILENAME=$(basename ${MITOGENOME/.fa/})

#now let's run mitofinder, the only option I added from the default is the --adjust-direction, which makes sure that all mitochondrias are outputted in the same direction
${MITOFINDER_PATH}/mitofinder -a ${MITOGENOME} \
	-t mitfi \
	-r ${REFERENCE} \
	-j ${FILENAME} \
	-o 2
