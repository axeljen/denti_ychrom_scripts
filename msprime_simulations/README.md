Scripts to run the msprime simulations as presented in the manuscript.

The directory autosomal_simulations contains python and wrapper scripts for running autosomal simulations, and the ychrom_simulations the same but for the ychromosomal simulations.

Executing the script autosomal_simulations/run_autosims.sh will simulate 100 chromosomes with the prespecified gene flow proportion from mitis to denti, and then calculate the Dstatistics for the simulated data using Dsuite.

Ychromosomal simulations are run by executing ychrom_simulations/run_ychrom_simulations.sh.

PAR simulation scripts and workflow are in par_simulations/.