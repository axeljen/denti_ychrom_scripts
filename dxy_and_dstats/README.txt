This directory contains the scripts and files used to calculate nucleotide divergence (dxy) and D-statistics.

Dxy was calculated between male samples among the mona and mitis groups, using the pixy_dxy.sh script together with the pops-denti-mitis-dxy-males.txt popfile.

D-stats was first calculated on a chrom-by-chrom basis with dsuite_dtrios.sh, and the chromosome-specific outputs were combined with dtrios_combine.txt. The file sets_samples.txt was used as the setsfile, and the all_samples_topology.nwk as the treefile.

The fd_windows.sh script was used to calculate Fd/dstatistics in 10k sliding windows along the genome, together with the fd_popfile.txt.