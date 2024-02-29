import msprime
import Bio
from Bio import bgzf
import sys

# take the proportion of gene flow from mitis to denti as a command line argument
gf_prop = float(sys.argv[1])
# out prefix as second
out_prefix = sys.argv[2]

ploidy = 1

############# These parameters are set manually, based on the demographic inference results #############

# generation time, for scaling times below to generations
g = 10
mu = 4.82e-9
# scale factor, for scaling the Ne's below when running on non-autosomal data
scale_factor = 4

# dictionary with Nes
Nes = {
	'denti': 200000 / scale_factor,
	'wolfi': 220000 / scale_factor,
	'pogonias': 230000 / scale_factor,
	'mona': 160000 / scale_factor,
	'neglectus': 110000 / scale_factor,
	'mitis': 510000 / scale_factor,
	'nictitans': 280000 / scale_factor,
	'cephus': 470000 / scale_factor,
	'denti_wolfi': 41000 / scale_factor,
	'denti_wolfi_pogonias': 250000 / scale_factor,
	'denti_wolfi_pogonias_mona': 200000 / scale_factor,
	'denti_wolfi_pogonias_mona_neglectus': 330000 / scale_factor,
	'denti_wolfi_pogonias_mona_neglectus_mitis_cephus_nictitans': 310000 / scale_factor,
	'mitis_nictitans': 130000 / scale_factor,
	'mitis_nictitans_cephus': 690000 / scale_factor,
	'macaque': 100000 / scale_factor,
	'root': 270000 / scale_factor,
}

# dictionary with split times
split_times = {
	'denti': 0,
	'wolfi': 0,
	'pogonias': 0,
	'mona': 0,
	'neglectus': 0,
	'mitis': 0,
	'nictitans': 0,
	'cephus': 0,
	'denti_wolfi': 1100000 / g,
	'denti_wolfi_pogonias': 1600000 / g,
	'denti_wolfi_pogonias_mona': 3300000 / g,
	'denti_wolfi_pogonias_mona_neglectus': 6200000 / g,
	'denti_wolfi_pogonias_mona_neglectus_mitis_cephus_nictitans': 7300000 / g,
	'mitis_nictitans': 2500000 / g,
	'mitis_nictitans_cephus': 3800000 / g,
	'root': 17000000 / g,
}

# dictionaries with gene flow events
gf = {
	1: {'dest': 'mitis', 'source': 'denti', 'time': 1000000 / g, 'proportion': gf_prop},
	2: {'dest': 'cephus', 'source': 'denti_wolfi_pogonias', 'time': 2700000 / g, 'proportion': 0.085},
	3: {'dest': 'mitis_nictitans', 'source': 'denti_wolfi_pogonias', 'time': 3100000 / g, 'proportion': 0.021},
}

# set the sequence length
seqlen = 10000
# set the recombination rate
r = 0
# set the number of replicates
replicates = 1000


def run_simulation(seqlen=10000, ploidy=2, r=1e-8, Nes={}, split_times={}, replicates=1, seed=None, gf={}):
	""" Run a simulation with the given parameters and return the tree sequence. """
	demography = msprime.Demography()
	# Add all the populations
	## Starting with the tips
	print("Adding extant populations")
	demography.add_population(name="denti", initial_size=Nes['denti'])
	demography.add_population(name="wolfi", initial_size=Nes['wolfi'])
	demography.add_population(name="pogonias", initial_size=Nes['pogonias'])
	demography.add_population(name="mona", initial_size=Nes['mona'])
	demography.add_population(name="neglectus", initial_size=Nes['neglectus'])
	demography.add_population(name="mitis", initial_size=Nes['mitis'])
	demography.add_population(name="nictitans", initial_size=Nes['nictitans'])
	demography.add_population(name="cephus", initial_size=Nes['cephus'])
	demography.add_population(name="macaque", initial_size=Nes['macaque'])
	# and then, one by one, the ancestral populations
	print("Adding ancestral populations")
	demography.add_population(name="denti_wolfi", initial_size=Nes['denti_wolfi'])
	demography.add_population(name="denti_wolfi_pogonias", initial_size=Nes['denti_wolfi_pogonias'])
	demography.add_population(name="denti_wolfi_pogonias_mona", initial_size=Nes['denti_wolfi_pogonias_mona'])
	demography.add_population(name="denti_wolfi_pogonias_mona_neglectus", initial_size=Nes['denti_wolfi_pogonias_mona_neglectus'])
	demography.add_population(name="denti_wolfi_pogonias_mona_neglectus_mitis_cephus_nictitans", initial_size=Nes['denti_wolfi_pogonias_mona_neglectus_mitis_cephus_nictitans'])
	demography.add_population(name="mitis_nictitans", initial_size=Nes['mitis_nictitans'])
	demography.add_population(name="mitis_nictitans_cephus", initial_size=Nes['mitis_nictitans_cephus'])
	demography.add_population(name="root", initial_size=Nes['root'])
	# add the events in chronological order
	## gene flow from mitis to denti
	demography.add_mass_migration(time=gf[1]['time'], source=gf[1]['source'], dest=gf[1]['dest'], proportion=gf[1]['proportion'])
	## split between denti and wolfi
	demography.add_population_split(time=split_times['denti_wolfi'], derived=['denti', 'wolfi'], ancestral='denti_wolfi')
	## split between denti/wolfi and pogonias
	demography.add_population_split(time=split_times['denti_wolfi_pogonias'], derived=['denti_wolfi', 'pogonias'], ancestral='denti_wolfi_pogonias')
	## split between mitis and nictitans
	demography.add_population_split(time=split_times['mitis_nictitans'], derived=['mitis', 'nictitans'], ancestral='mitis_nictitans')
	## gene flow from cephus ancestor to denti/wolfi/pogonias ancestor
	demography.add_mass_migration(time=gf[2]['time'], source=gf[2]['source'], dest=gf[2]['dest'], proportion=gf[2]['proportion'])
	## gene flow from mitis/nictitans ancestor to denti/wolfi/pogonias ancestor
	demography.add_mass_migration(time=gf[3]['time'], source=gf[3]['source'], dest=gf[3]['dest'], proportion=gf[3]['proportion'])
	## split between denti/wolfi/pogonias and mona
	demography.add_population_split(time=split_times['denti_wolfi_pogonias_mona'], derived=['denti_wolfi_pogonias', 'mona'], ancestral='denti_wolfi_pogonias_mona')
	## split between cephus and mitis/nictitans
	demography.add_population_split(time=split_times['mitis_nictitans_cephus'], derived=['cephus', 'mitis_nictitans'], ancestral='mitis_nictitans_cephus')
	## split between denti/wolfi/pogonias/mona and neglectus
	demography.add_population_split(time=split_times['denti_wolfi_pogonias_mona_neglectus'], derived=['denti_wolfi_pogonias_mona', 'neglectus'], ancestral='denti_wolfi_pogonias_mona_neglectus')
	## split between denti/wolfi/pogonias/mona/neglectus and mitis/nictitans/cephus
	demography.add_population_split(time=split_times['denti_wolfi_pogonias_mona_neglectus_mitis_cephus_nictitans'], derived=['denti_wolfi_pogonias_mona_neglectus', 'mitis_nictitans_cephus'], ancestral='denti_wolfi_pogonias_mona_neglectus_mitis_cephus_nictitans')
	## split between denti/wolfi/pogonias/mona/neglectus/mitis/nictitans/cephus and macaque
	demography.add_population_split(time=split_times['root'], derived=['denti_wolfi_pogonias_mona_neglectus_mitis_cephus_nictitans', 'macaque'], ancestral='root')
	# set recombination rate
	# run the actual simulation
	ts = msprime.sim_ancestry(
		recombination_rate=r,
		sequence_length=seqlen,  
		samples={'denti':1,'wolfi':1,'pogonias':1,'mona':1,'neglectus':1,'mitis':1,'nictitans':1,'cephus':1,'macaque':1},
		demography = demography,
		random_seed=seed,
        ploidy=ploidy,
		num_replicates=replicates,
	)
	return ts

# run the simulation
ts = run_simulation(seqlen=seqlen, ploidy=ploidy, r=r, Nes=Nes, split_times=split_times, replicates=replicates, seed=None, gf=gf)


labels = {
	0: 'denti',
	1: 'wolfi',
	2: 'pogonias',
	3: 'mona',
	4: 'neglectus',
	5: 'mitis',
	6: 'nictitans',
	7: 'cephus',
	8: 'macaque',
}

# list to store the trees
monophyletic_trees = []
# dictionary to keep track of distances
distances = {}

# iterate over the trees
for i,t in enumerate(ts):
	for tree in t.trees():
		# get the distance between denti and mitis
		if ploidy == 1:
			dist = tree.path_length(0,5)
		#elif ploidy == 2:
		#	dist = tree.path_length(0,10)
		# add the distance to the dictionary
		if dist in distances:
			distances[dist] += 1
		else:
			distances[dist] = 1
		# if the distance is more than 3, add the tree to the list
		if ploidy == 1 and dist == 2:
			monophyletic_trees.append(tree.newick(include_branch_lengths=False,
							node_labels=labels,))
		elif ploidy == 2 and dist < 5:
			monophyletic_trees.append(tree.newick(include_branch_lengths=False,
							node_labels=labels,))

# sort the distance dictionary by distance
distances = dict(sorted(distances.items()))

# write the counts and trees to file
with open(out_prefix + "_distance_counts.txt", "w") as f:
	f.write("Distance\tCount\n")
	for dist in distances:
		f.write(str(dist) + "\t" + str(distances[dist]) + "\n")
with open(out_prefix + "_trees.txt", "w") as f:
	f.write("\n".join(monophyletic_trees))


