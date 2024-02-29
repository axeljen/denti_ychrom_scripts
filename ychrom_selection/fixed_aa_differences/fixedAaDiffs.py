import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.append('/domus/h1/axeljen/phylogenomics')
import functions as fn
import argparse

# define a function to convert a fasta sequence to an amino acid sequence
def translate(seq, break_on_stop = False):
	# make a list of codons by splitting the sequence in triplets
	codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
	aaseq = ''
	# loop over the codons and translate them
	for codon in codons:
		aa = fn.codon2aa(codon)
		if aa == "*" and break_on_stop:
			break
		aaseq += aa
	return aaseq


# function to check if the two groups are fixed for a change
def checkFixedChange(pops, aln, pos):
	# fetch all alleles at this position for both groups
	alleles = {}
	missingness = {p: 0 for p in pops}
	denti_missing = False
	mitis_missing = 0
	for pop in pops:
		alleles[pop] = []
		for sample in pops[pop]:
			allele = aln.sequences[sample].sequence[pos]
			if not allele in ['-', '?', 'X']:
				alleles[pop].append(allele)
			else:
				missingness[pop] += 1
				if sample == 'FK104_C_denti':
					denti_missing = True
				elif pop == "blue":
					mitis_missing += 1
		# if either pop only has missing data, return false
		if missingness[pop] == len(pops[pop]):
			return False,'missing', alleles
		# same if either denti or all mitis are missing
		if denti_missing or mitis_missing == len(pops['blue']) - 1:
			return False, 'missing', alleles
	# check if the two groups are fixed for a change
	for allele in alleles[list(pops.keys())[0]]:
		# if the allele is in the other group, return false
		if allele in alleles[list(pops.keys())[1]]:
			return False, 'shared', alleles
	# if we get here, the two groups are fixed for a change
	return True, 'fixed', alleles


# prep a parser for the command line arguments
parser = argparse.ArgumentParser(description='Identify/count fixed amino acid changes between two groups, based on either an aa or nt sequence.')

# popfile with samples and groups for the samples to consider
parser.add_argument('-p', '--popfile', type=str, required=True, help='Popfile with samples and groups for the samples to consider.')
# fasta file with sequences
parser.add_argument('-a', '--alignment', type=str, required=True, help='Alignment file with sequences.')
# output file
parser.add_argument('-o', '--output', type=str, required=True, help='Output file.')
# store true argument to check if we should translate or not
parser.add_argument('-t', '--translate', action='store_true', help='Translate the sequences to amino acids before comparing them.')
# path to output alignment file if we wanna write it
parser.add_argument('-w', '--write-aln', type=str, help='Path to output alignment file if we wanna write it.', default = None)

# parse the arguments
args = parser.parse_args()

# read the popfile
pops = fn.parsePopfile(args.popfile)
# and make a list of all the samples to keep from this
samples = []
for pop in pops:
	samples += pops[pop]

# read the alignment 
aln = fn.readSequenceFile(args.alignment)

# prune out the samples that are not in the popfile
aln.subsetSamples(samples)

# if translate is true, do this now
if args.translate:
	aln.nt2aa()

# make a dictionary to store the fixed changes
fixed_changes = {}
fixed_changes_count = 0
missing_data = 0
evaluated_sites = 0

# loop through all the positions in the alignment
for pos in range(aln.length):
	# check if the two groups are fixed for a change
	fixed, reason, alleles = checkFixedChange(pops, aln, pos)
	# if they are, add this to the dictionary
	if fixed:
		fixed_changes[pos] = (reason, alleles)
		fixed_changes_count += 1
		evaluated_sites += 1
	elif reason == 'missing':
		missing_data += 1
	elif reason == 'shared':
		evaluated_sites += 1

# print the results
print('Evaluated sites: %i' % evaluated_sites)
print('Fixed changes: %i' % fixed_changes_count)
print('Missing data: %i' % missing_data)

# write the output
with open(args.output + '_fixed_alleles.txt', 'w') as outfile:
	outfile.write('pos\t{}_allele\t{}_allele\n'.format(list(pops.keys())[0], list(pops.keys())[1]))
	for site in fixed_changes:
		a1 = '/'.join(set(fixed_changes[site][1][list(pops.keys())[0]]))
		a2 = '/'.join(set(fixed_changes[site][1][list(pops.keys())[1]]))
		outfile.write('{}\t{}\t{}\n'.format(site, a1, a2))

# write the counts, too
with open(args.output + '_fixed_diffs_counts.txt', 'w') as outfile:
	outfile.write('good_sites\tfixed_changes\tmissing_data\n')
	outfile.write('{}\t{}\t{}\n'.format(evaluated_sites, fixed_changes_count, missing_data))

# write the output alignments
if args.write_aln:
	fn.writeSequenceFile(aln, args.write_aln)

	


