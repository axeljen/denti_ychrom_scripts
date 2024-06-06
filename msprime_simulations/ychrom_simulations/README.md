# Y-chromosomal simulations

This directory contains scripts to simulate a non-recombining, haploid loci mimicking the Y chromosome. Simulations are run inside the denti_simulations_ychrom.py script, which contains all simulation parameters etc, just takes a migration rate from mitis to denti as an argument.

The scripts also counts the number of times denti and mitis are monophyletic.

The file array_proportions.txt just contains the geneflow proportions (in second column) to use for each array in the job (first column is the array number).

The script run_ychrom_simulations.sh runs simulations in an slurm array job, one task for each migration rate, each task running one hundred replicate simulations, and output the monophyly counts to a file.

## Dependencies
These scripts use the following software and packages (versions used in parentheses):
- msprime (1.2.0)
- tidyverse (2.0.0)

The simulations are simply run as follows:

```bash
# the script is set to output files to the out/ directory, so make sure it exists
mkdir -p out

# give an output name as argument
sbatch run_ychrom_simulations.sh msprime_ychrom_simulations
```

When the simulations are done, the plotYsims.R script can be used to plot the monophyly counts (as shown in Figure S11)

```bash
Rscript plotYsims.R
```
