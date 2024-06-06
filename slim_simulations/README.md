Scripts to run the forward simulations in SLIM. The model is specified in simple_yintro_model.slim. 

The simulations are run in the slimsim.sh, and automatated slurm submissions for different migration rates etc. are controlled through the slim_wrapper.sh script.

## Dependencies

These scripts use the following software and packages (versions used in parentheses):
- SLIM (4.0.1)
- tidyverse (2.0.0)

## Workflow

1. Run the sumulations by executing slim_wrapper.sh

```bash
bash slim_wrapper.sh
```

2. When all simulations are finished, concatenate all the output files to a single file:

```bash
# name of output directory
out=simple_yintro_model

# start with the header of the first file
head -n 1 $(find ${out} -name "*.txt" | head -n 1) > ${out}/simulations_merged.txt
# append the rest of the files
for f in $(find ${out} -name "*txt" | grep -v simulations_merged.txt); do
	tail -n +2 $f >> ${out}/simulations_merged.txt
done

```

3. Plot output with plotSlimSims.R, which will recreate Figure 6C.

```bash

Rscript plotSlimSims.R

```