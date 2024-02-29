Scripts used to run the bpp-msci models.

Controlfiles for the two models included in the manuscript are included, where model 1 is presented in the main paper and 2 in the supplementary material.

I'd run two independent runs for each model, and then check convergence and merge them.

the run_bpp.sh script is just the execution script to actually run bpp.

When the runs finished, I used bpp_parse.R to combine, scale and plot the output files. Note that this script is written for this particular usecase, changing the number of hybridization nodes etc is very likely to break it.