# wrapper around the sbatch submission of slim simulations

# popsize to simulate
POPSIZE=200000

# name of simulation
NAME=simple_yintro_model

# slim script to run
#SLIMSCRIPT=simple_yintro_model_v1.slim
SLIMSCRIPT=simple_yintro_model.slim

# outdir
OUTDIR=${NAME}

# loop over replicates
for i in $(seq 1 ${NREPS})
do
# total introproportions/inital allele frequencies (Y-chromosomal frequency will be twice this)
for INTROPROP in 0.0005 0.001 0.002 0.003 0.004 0.005
do
# for each allele frequency, loop over selection coefficients
for SEL in 0 0.00001 0.0001 0.001 0.01
do
echo ${INTROPROP} ${SEL}
# submit job
sbatch slimsim.sh ${SLIMSCRIPT} ${SEL} ${POPSIZE} ${INTROPROP} ${NAME}
done
done
done