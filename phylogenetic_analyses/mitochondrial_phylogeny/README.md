
This directory contains the scripts and input files used to construct the mitochondrial phylogeny.

## Dependencies

These scripts use the following software and packages (versions used in parentheses):
- iqtree (2.2.2.6 )

## Workflow

The mitophylogeny.sh script take the alignment and coordinates of partitions to 1. perform modeltest to find best model for each partition and 2. contruct a phylogeny using these models.

```bash
# submit script with output name as argument
sbatch mitophylogeny.sh mitophylogeny
```