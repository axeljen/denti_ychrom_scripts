
Scripts used to run the branch site model with in codeml-paml.

All the input parameters are specified in the codeml_wrapper.sh script, which modifies the codeml_slurm_raw.sh script into a new one, called codeml_slurm.sh.

This script then runs the codeml.py script which is a python wrapper to codeml.

## Dependencies

These scripts use the following software and packages (versions used in parentheses):
- PAML (4.9j)
- my phylogenomics repo: https://github.com/axeljen/phylogenomics. After cloning and installing dependencies, the path to this repo should be modified in the codeml.py script.

## Workflow

Specify the parameters inside the codeml_wrapper.sh script and run as:

```bash
bash codeml_wrapper.sh
```
