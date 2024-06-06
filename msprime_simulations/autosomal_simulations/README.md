# Autosomal simulations

This directory contains scripts to run autosomal simulations in msprime, and estimate D-statistics on these simulations, to recreate Figure 6A.

The actual msprime simulation parameters are specified in the denti_simulations_autosomes.py script. The Dstat_SETS.txt and Dstat_tree.txt contains population assignments and topology for the D-statistic calculations, respectively. parse_dsuite_output.py just parses the output from Dsuite and makes a big table of everything for plotting.

## Dependencies
These scripts use the following software and packages (versions used in parentheses):
- msprime (1.2.0)
- Dsuite (0.4 r42)
- tidyverse (2.0.0)

## Workflow

1. Loop over migration rates from mitis to denti and execute the simulation script run_autosomes.sh, which takes the migration rate as an argument.

```bash

# output file name
outname=autosomal_sims

# number of replicates to run
NREPS=10

for i in $(seq 1 ${NREPS}); do
for mig in $(seq 0 0.0005 0.01); do
echo ${i} ${mig}
sbatch run_autosims.sh ${mig} ${outname}
done
done

# when these are done, plot the output with plot_simulated_D.R (make sure that the input file is adjusted correctly inside script), to recreate Figure 6A.
Rscript plot_simulated_D.R

```
