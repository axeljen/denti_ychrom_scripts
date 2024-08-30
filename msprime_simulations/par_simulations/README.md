# PAR simulations

This directory contains the scripts used to simulate a Y-chrom with a linked pseudoautosomal region (PAR) and autosomes under varying PAR recombination rate, to explore the expected introgression on PAR in different scenarios.

## Dependencies

These scripts require the following dependencies:
-msprime
-Simon Martin's github repo (https://github.com/simonhmartin/genomics_general)

## Workflow

The following workflow is used to simulate Y-chrom introgression with PAR and autosomes under the same, tenfold, and twentyfold recombination rate on PAR relative to autosomes:

```bash

# loop through the recombination rates
for rate in 1 10 20
do
# runsims.sh will run the simulate_parintro.py script and submit a script for estimating fd
sbatch runsims.sh $rate
done


```
