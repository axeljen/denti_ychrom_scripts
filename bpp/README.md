Scripts used to run the bpp-msci models.

Controlfiles for the two models included in the manuscript are included, where model 1 is presented in the main paper and 2 in the supplementary material.

I'd run two independent runs for each model, and then check convergence and merge them.

the run_bpp.sh script is just the execution script to actually run bpp.

When the runs finished, I used bpp_parse.R to combine, scale and plot the output files. Note that this script is written for this particular usecase, changing the number of hybridization nodes etc is very likely to break it.

## Dependencies

For these scripts to work, BPP needs to be installed and in PATH. The first steps are also using my python parsing scripts, which can be cloned from here:
https://github.com/axeljen/phylogenomics.

In addition to Python3 (https://www.python.org/downloads/), pysam and biopython are required for these scripts:

After installing python3, install the required packages with pip:
	python3 -m pip install pysam biopython


## Workflow description

1. Sample intergenic loci from the reference genome and extract to alignments.

I did this as described in sample_loci.sh. Note that this is a rather bulky script just documenting the steps I took interactively on the command line, so they need to be ran step-by-step. The output of this step is the 1000_loci_bpp_format.txt file, containing the alignments in the format required by bpp.

2. Create the control file to use with bpp. I followed the bpp documentation and created two msci models, each with three migration bands but model_2 has the timing of gene-flow event 1 and 2 switched around. The control files I used are included in this repository.

3. Run bpp. I ran two independent runs per model on a slurm cluster:
```
	mkdir model_1_run1 model_1_run2 model_2_run1 model_2_run2
	cp model_1.ctl model_1_run1
	sbatch run_bpp.sh model_1_run1/model_1.ctl
	cp model_1.ctl model_1_run2
	sbatch run_bpp.sh model_1_run2/model_1.ctl
	cp model_2.ctl model_2_run1
	sbatch run_bpp.sh model_2_run1/model_2.ctl
	cp model_2.ctl model_2_run2
	sbatch run_bpp.sh model_2_run2/model_2.ctl
```
4. bpp_parse.R was used to compare, join and summarize the mcmc output from bpp, and to plot the results. Note that this script is very specific to this usecase, and visualizing the output of bpp can be done in many ways. The script needs the following packages (versions I used are indicated here): tidyverse/2.0.0, bppr/0.6.3, ape/5.8, coda/0.19-4.1, ggpubr/0.6.0, psych/2.4.3, MoreTreeTools/0.0.1, treeio/1.26.0, ggrepel/0.9.5.


