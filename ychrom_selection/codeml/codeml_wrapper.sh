# wrapper around the python-codeml script

# name of the test
name=bsm_denti-mitis


# path to alignmnent to analyze
aln=../gene_alignments/KDM5D_XP_014984073.1.phy
# tree
tree="tree_allmales.nwk"
# samples 
samples="samples.txt"
# outgroup
outgroup="SAMN02692280_Chlorocebus_sabaeus_saintkitts"
# foreground branches
foreground="FK104_C_denti,PD_0096_Cercopithecus_mitis"

# we'll need to sed into the codeml raw to set the foreground branches
cat codeml_slurm_raw.sh | sed "s/FOREGROUND_BRANCHES/${foreground}/g" > codeml_slurm.sh

# submit to slurm
sbatch codeml_slurm.sh ${aln} ${tree} ${samples} ${name} ${outgroup}
