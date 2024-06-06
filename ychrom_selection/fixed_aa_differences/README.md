
This directory contains the scripts used to identify fixed amino acid differences between the denti/mitis branch and the wolfi/pogonias branch on the Y chromosome.

## Dependencies
- my phylogenomics repo: https://github.com/axeljen/phylogenomics. After cloning and installing dependencies, the path to this repo should be modified in the fixedAaDiffs.py script.

## Workflow

- Using a population assignments file (two columns, sample\tpop) the script fixedAaDiffs.py can be used to find fixed AA-substitutions between these populations/groups using an inframe, transcript alignment.

After setting all parameters in the find_diffs_wrapper.sh, simply execute the wrapper script:

```bash
bash find_diffs_wrapper.sh
```

This will create a file called name.counts_incl_stops.txt, which is equivalent to Table S5.