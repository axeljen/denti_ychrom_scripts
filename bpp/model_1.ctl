seed =  -1

seqfile = 1000_loci_bpp_format.txt
Imapfile = msci_imap.txt
outfile = msci_model_1.txt
mcmcfile = msci_model_1.mcmc

# fixed number of species/populations 
speciesdelimitation = 0

# fixed species tree
speciestree = 0

species&tree = 9  C.denti  C.wolfi  C.pogonias  C.mona  C.neglectus  C.cephus  C.nictitans  C.mitis  M.mulatta
                  1  1  1  1  1  1  1  1  1
                 (((((((((C.denti)B[&phi=0.800000,tau-parent=yes],C.wolfi)denti_wolfi,C.pogonias)Mon1)Q[&phi=0.800000,tau-parent=yes])W[&phi=0.800000,tau-parent=yes],C.mona)Mon2,C.neglectus),((Q[&phi=0.200000,tau-parent=no],C.cephus)P,(W[&phi=0.200000,tau-parent=no],(C.nictitans,(B[&phi=0.200000,tau-parent=no],C.mitis)A)mitis_anc)Z)mitis_cephus_anc)guenon_root, M.mulatta);

# unphased data for all guenons, haploid sequence for macaque
phase =   1  1  1  1  1  1  1  1  1

# use sequence likelihood
usedata = 1

nloci = 1000

Threads = 20 1 1

# do not remove sites with ambiguity data
cleandata = 0

thetaprior = 3 0.006 e # gamma(a, b) for theta (estimate theta)
tauprior = 3 0.04 # gamma(a, b) for root tau & Dirichlet(a) for other tau's
phiprior = 1 1

# finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
finetune =  1: 0.01 0.02 0.03 0.04 0.05 0.01 0.01

# MCMC samples, locusrate, heredityscalars, Genetrees
print = 1 0 0 0   * 
burnin = 20000
sampfreq = 2
nsample = 200000
