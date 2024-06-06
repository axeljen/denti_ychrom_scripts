#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 1-00:00:00
#SBATCH -J fd_windows
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# load modules
module load bioinfo-tools bcftools

# give path to Simon Martin's scripts (https://github.com/simonhmartin/genomics_general)
SCRIPT_PATH=/home/axeljen/genomics_general

# give the full path to vcf as command line input
INPUT_VCF=${1}

# will write to ./out, make sure it exists
mkdir -p out

# popfile as second 
POPS=fd_popfile.txt

# specify taxa to test
P1=wolfi
P2=denti
P3=mitis
OUT=mmul

# window size and step
WINDOW=10000
STEP=10000

# How many cpu threads are we running on?
THREADS=16

# will extract chromomsome info with bcftools, assuming this is a single chrom vcf
CHROM=$(bcftools view -H ${INPUT_VCF} | head -n 1 | cut -f 1)

# then we'll convert this to a geno file, put it on tmpstorage
GENO=${SNIC_TMP}/${SLURM_JOB_ID}_${CHROM}.geno.gz
# use Simon's parseVCF.py for this
python3 ${SCRIPT_PATH}/VCF_processing/parseVCF.py -o ${GENO} -i ${INPUT_VCF}

# let's create an output path from the info above
OUTFILE=out/${CHROM}_p1-${P1}_p2-${P2}_p3-${P3}_out-${OUT}_w-${WINDOW}_step-${STEP}.csv.gz
# and cat the parameters, including the popfile, to a file next to this one so that we keep this info
echo -e "####\n\nP1: ${P1}; P2: ${P2}; P3: ${P3}; Outgroup: ${OUT};\nWindow size: ${WINDOW};\nStep size: ${STEP};\n\n####\nPopfile used:\n\n" > ${OUTFILE/.csv.gz/params.txt}
cat ${POPS} >> ${OUTFILE/.csv.gz/params.txt}

# then run the command
python3 ${SCRIPT_PATH}/ABBABABAwindows.py -g ${GENO} -o ${OUTFILE} -f phased \
	-P1 ${P1} -P2 ${P2} -P3 ${P3} -O ${OUT} \
	--windType coordinate \
	-w ${WINDOW} -s ${STEP} \
	--popsFile ${POPS} \
	--writeFailedWindows \
	-T ${THREADS}





