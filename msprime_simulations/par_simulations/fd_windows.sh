#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p node
#SBATCH -N 20
#SBATCH -t 1-00:00:00
#SBATCH -J fd_windows
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# give path to Simon Martin's scripts
SCRIPT_PATH=/home/axeljen/genomics_general/

# give the full path to vcf as command line input
INPUT_VCF=$1

# rho factor on par (only used to include this in the output table for simple plotting)
PAR_RHO_FACTOR=$2

# gene flow proportion (only used to include this in the output table for simple plotting)
gf_prop=$3

# popfile as second 
POPS=sets.txt

# P1-outgroup next
P1=wolfi
P2=denti
P3=mitis
OUT=macaca

# window size and step
WINDOW=50000
STEP=50000

# threads
THREADS=20

# then we'll convert the vcf file to geno, put it on tmpstorage
GENO=${SNIC_TMP}/${SLURM_JOB_ID}.geno.gz
# use Simon's parseVCF.py for this
python3 ${SCRIPT_PATH}/VCF_processing/parseVCF.py -o ${GENO} -i ${INPUT_VCF}

# make output file variable
OUTFILE=out/$(basename ${INPUT_VCF%%.vcf.gz}).csv.gz

# then run the command
python3 ${SCRIPT_PATH}/ABBABABAwindows.py -g ${GENO} -o ${OUTFILE} -f phased \
	-P1 ${P1} -P2 ${P2} -P3 ${P3} -O ${OUT} \
	--windType coordinate \
	-w ${WINDOW} -s ${STEP} \
	--popsFile ${POPS} \
	--writeFailedWindows \
	-T ${THREADS}

# when done, add the rho factor as an additional column for convenience
zcat ${OUTFILE} | head -n 1 | awk ' {print $0",rho_factor,gf_prop"} ' > tmp_$(basename ${OUTFILE%%.gz})
zcat ${OUTFILE} | tail -n +2 | awk -v rho=${PAR_RHO_FACTOR} -v gf_prop=${gf_prop} ' {print $0","rho","gf_prop} ' >> tmp_$(basename ${OUTFILE%%.gz})

mv tmp_$(basename ${OUTFILE%%.gz}) ${OUTFILE%%.gz}


