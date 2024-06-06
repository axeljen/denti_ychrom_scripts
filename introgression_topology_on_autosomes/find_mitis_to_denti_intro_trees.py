from ete3 import Tree
import os
import argparse
import sys
import time

def root_tree(tree,outgroup,verbose=False):
	try:
		if len(outgroup) > 1:
			out_ancestor = tree.get_common_ancestor(outgroup)
			tree.set_outgroup(out_ancestor)
		else:
			tree.set_outgroup(outgroup[0])
	except:
		if verbose:
			print("Could not root tree.")
		tree = "NA"
	return tree

def check_monophyly(tree, samples):
	if tree.check_monophyly(values = samples, target_attr="name")[0]:
		return True
	else:
		return False

# parse arguments
parser = argparse.ArgumentParser(description='Find trees with introgression from mitis to denti')
parser.add_argument('-t', '--treefile', help='File with trees', required=True)
parser.add_argument('-s', '--samples', help='File with samples and populations', required=True)
parser.add_argument('-o', '--output', help='Output file', required=True)
args = parser.parse_args()

treefile = args.treefile

samples = args.samples

output = args.output

# parse popfile
popdict = {}
with open(samples) as f:
	for line in f.readlines():
		if line.strip() == "":
			continue
		sample, pop = line.strip().split()
		popdict[sample] = pop

# outgroup taxon for rooting trees
outgroup = "macaca"

trees = []

with open(treefile) as f:
	for line in f.readlines():
		if line.strip().startswith("window"):
			# take the headers and then skip line
			headers = line.strip().split()
			continue
		values = line.strip().split()
		ldict = dict(zip(headers,values))
		if ldict.get('tree') is None:
			continue
		if ldict['tree'] == "NA":
			continue
		else:
			# parse the tree
			t = Tree(ldict.get('tree'))
			if t is None:
				print("No tree found.")
				continue
			# prune
			t.prune([s for s in popdict.keys()])
			# root on outgroup if possible
			tree = root_tree(t, [o for o in popdict.keys() if popdict[o] == outgroup])
			# next, check if denti is monophyletic
			if not check_monophyly(tree, [s for s in popdict.keys() if popdict[s] == "denti"]):
				continue
			# check if denti and mitis is monophyletic
			if not check_monophyly(tree, [s for s in popdict.keys() if popdict[s] == "denti" or popdict[s] == "mitis"]):
				continue
			# check if mitis/denti/nicitans is monophyletic
			if not check_monophyly(tree, [s for s in popdict.keys() if popdict[s] == "denti" or popdict[s] == "mitis" or popdict[s] == "nictitans"]):
				continue
			# check if mona, pogonias and wolfi are monophyletic
			if not check_monophyly(tree, [s for s in popdict.keys() if popdict[s] == "mona" or popdict[s] == "pogonias" or popdict[s] == "wolfi"]):
				continue
			trees.append({
				'chrom': ldict['chrom'],
				'start': ldict['start'],
				'end': ldict['end'],
				'sites':ldict['sites'],
				'tree': tree.write(format=5),
				'tree_full': ldict['tree'],
			})

if len(trees) == 0:
	print("No trees found")
with open(output, 'w') as f:
	f.write("chrom\tstart\tend\tsites\ttree\ttree_full\n")
	for tree in trees:
		f.write(f"{tree['chrom']}\t{tree['start']}\t{tree['end']}\t{tree['sites']}\t{tree['tree']}\t{tree['tree_full']}\n")