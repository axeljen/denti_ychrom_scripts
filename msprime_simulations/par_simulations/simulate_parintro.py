import msprime
import sys
import pysam

## script to simulate introgression on a chromosome where one part is non-recombining, to see how the 
# recombining part is affected by gene flow

# prefix of output files as first argument
outprefix = sys.argv[1]

# number of repetitions as second argument
#nreps = 100
nreps = int(sys.argv[2])

# rho factor on par as command line input
#rho_factor = 10
rho_factor = float(sys.argv[3])

# gene flow proportion from dmitis to denti
#gf_prop = 0.01
gf_prop = float(sys.argv[4])
# initiate a demography event 
guenons = msprime.Demography()

# using a single Ne across both species
Ne = 50000

# recombination rate on the recombining part "autosome" and par, using the estimated macaque rate of 4.48e-9
recomb_rate = 4.48e-9

# sequence length to sequence
seq_len = 350000

# ratemap where the first 50000 bp has no recombination (Ychrom), then a par with "regular" recombination rate, than three "autosomes" separated by 1 bp with recombination rates of 0.5.
ratemap = msprime.RateMap(position = [0,50000,200000,200001,250000,250001,300000,300001,350000],rate = [0,recomb_rate * rho_factor,0.5,recomb_rate,0.5,recomb_rate,0.5,recomb_rate])

# mutation rate
mu=4.82e-9

# generation time, for converting time to generations
g = 10


# add tips
guenons.add_population(name="denti", initial_size=Ne)
guenons.add_population(name="wolfi", initial_size=Ne)
guenons.add_population(name="mitis", initial_size=Ne)
guenons.add_population(name="macaca", initial_size=Ne)

# add internal nodes, sorted by age
guenons.add_population(name="denti_wolfi", initial_size=Ne)
guenons.add_population(name="denti_wolfi_mitis", initial_size=Ne)
guenons.add_population(name="denti_wolfi_mitis_macaca", initial_size=Ne)

# gene flow between denti and mitis
guenons.add_mass_migration(time = 1e6 / g, source = "denti", dest = "mitis", proportion = gf_prop)

# add splits, in millions of years
guenons.add_population_split(time = 2e6 / g, derived=["denti", "wolfi"], ancestral="denti_wolfi")
guenons.add_population_split(time = 7e6 / g, derived=["denti_wolfi", "mitis"], ancestral="denti_wolfi_mitis")
guenons.add_population_split(time = 12e6 / g, derived=["denti_wolfi_mitis", "macaca"], ancestral="denti_wolfi_mitis_macaca")


# dictionary with sample labels
labels = {
	0: 'denti',
	2: 'wolfi',
	3: 'mitis',
	4: 'macaca',
	}

# simulate 100 autosomes
sims_auto = msprime.sim_ancestry(
		recombination_rate=ratemap,
		sequence_length=seq_len,  
		samples={'denti': 1, 'wolfi': 1, 'mitis': 1, 'macaca': 1},
		demography = guenons,
		random_seed=None,
		ploidy=2,
		num_replicates = nreps,
		record_migrations = True,
	)

# function to extract migration tracts
def get_migration_tracts(ts):
	mitis_id = [p.id for p in ts.populations() if p.metadata['name'] == 'mitis'][0]
	migrating_tracts = []
	# get Y-tracts going directly from denti to mitis
	for m in ts.migrations():
		if m.dest == mitis_id:
			if m.left <= 50000:
				migrating_tracts.append((m.left,m.right))
	return migrating_tracts

# simulate mutations
for i,ts in enumerate(sims_auto):
	# get any remaining migration tracts from mits in denti
	mt = get_migration_tracts(ts)
	# if the Y remains in the population, sim mutations and write vcf output
	if len(mt) > 0:
		muts = msprime.sim_mutations(ts, rate = mu)
		with open(outprefix + "_longpar_auto_{}_{}.vcf".format(i, gf_prop), "w") as vcf_file:
			muts.write_vcf(vcf_file, contig_id="sim_{}_{}".format(i, rho_factor),
			individual_names = ["denti","wolfi","mitis","macaca"],
			)
	sys.stdout.write("Finished with simulation {}.".format(i))
