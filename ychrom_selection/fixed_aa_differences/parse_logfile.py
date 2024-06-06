import sys

name = sys.argv[1] #

logfile = name + ".log"
counts = name + ".all_counts.txt"

# fetch genes and counts
genes = []

for line in open(counts):
	if line.startswith("gene"):
		header = line.strip().split("\t")
		continue
	else:
		values = line.strip().split("\t")
		d = dict(zip(header, values))
		genes.append(d)


current_gene = None
current_sample = None
with open(logfile) as f:
	for line in f:
		if line.strip() == "":
			continue
		if line.strip() in [i['gene'] for i in genes]:
			if current_gene is None:
				stop_codons = []
				current_gene = line.strip()
			else:
				[d for d in genes if d['gene'] == current_gene][0]['stop_codons'] = stop_codons
				stop_codons = []
				current_gene = line.strip()
		if 'Translating' in line.strip():
			current_sample = line.strip().split(" ")[2]
		if 'stop' in line.strip(" "):
			if not current_sample in stop_codons:
				stop_codons.append(current_sample)
	# add last gene
	[d for d in genes if d['gene'] == current_gene][0]['stop_codons'] = stop_codons

for gene in genes:
	try:
		print(gene['gene'], gene['stop_codons'])
	except KeyError:
		print(gene['gene'], "no stop codons")
		
# write to file
with open(name + ".counts_incl_stops.txt", "w") as f:
	f.write('gene\tgoodsites\tfixed_diffs\tmissing_sites\tstop_codons\n')
	for gene in genes:
		f.write(gene['gene'] + "\t" + gene['goodsites'] + "\t" + gene['fixed_diffs'] + "\t" + gene['missing_sites'] + "\t" + ';'.join(gene['stop_codons']) + "\n")
