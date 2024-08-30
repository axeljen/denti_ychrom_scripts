#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 0-10:00:00
#SBATCH -J parsims
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

module load bioinfo-tools bcftools

# make sure the outdir and vcf dir exists
mkdir -p out
mkdir -p vcf_files

# factor increase in recombination rate on par as command line input
par_rho_factor=$1

# gene flow proportion
gf_prop=0.01

# number of replicates to run (these will be output as separate 'genomes'/'chromosomes' but in the same vcf file)
nreps=100

# output vcf path
out_vcf=vcf_files/rho_factor_${par_rho_factor}_gf_${gf_prop}.vcf.gz

# submit the fd script in a separate job with dependency
sbatch --dependency=afterok:${SLURM_JOB_ID} fd_windows.sh ${out_vcf} ${par_rho_factor} ${gf_prop}

# make a temporary directory to store the vcf files
mkdir tmp_${SLURM_JOB_ID}

# run simulations and put vcf files in the temporary directory
python3 simulate_parintro.py tmp_${SLURM_JOB_ID}/par_rho_factor_${par_rho_factor} ${nreps} $par_rho_factor $gf_prop

# when this is done, concatenate the vcf files into one
bcftools concat -Oz -o ${out_vcf} $(find tmp_${SLURM_JOB_ID} -name "*.vcf")

# remove the tempdir
rm -r tmp_${SLURM_JOB_ID}
