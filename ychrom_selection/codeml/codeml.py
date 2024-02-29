# import codeml and chi2 from bio
from Bio.Phylo.PAML import codeml
from Bio.Phylo.PAML.chi2 import cdf_chi2
from ete3 import Tree
import sys
sys.path.append("/home/axeljen/phylogenomics/")
import functions as fn
import re
import importlib
import os
import argparse

class testargs:
	def __init__(self, alignment, tree, samples, outgroup, foreground, maxmissing, allow_stop_codons, workdir, name, df):
		self.alignment = alignment
		self.tree = tree
		self.samples = samples
		self.outgroup = outgroup
		self.foreground = foreground
		self.maxmissing = maxmissing
		self.allow_stop_codons = allow_stop_codons
		self.workdir = workdir
		self.name = name
		self.df = df

def prep_tree(tree, foreground_clades, outgroup=None):
	# set root if given
	if outgroup is not None:
		tree.set_outgroup(outgroup)
	# and now we want to remove the branchlengths around the root 
	root = tree.get_tree_root()
	for child in root.get_children():
		if not child.is_leaf():
			child.delete()
	# add temporary labels to these nodes
	for node in tree.traverse():
		# get all the leaves of this node
		if not node.is_leaf():
			leaves = node.get_leaves()
			for clade in foreground_clades:
				if set(clade) == set([i.name for i in leaves]):
					node.add_feature("label","foreground")
		else:
			if [node.name] in foreground_clades:
				node.add_feature("label","foreground")
	# make the treestring
	treestring = re.sub("\[&&NHX:label=foreground\]", " #1",tree.write(format=9, features=["label"]))
	print(treestring)
	return treestring

# read alignment file
def prep_alignment(aln, samples, maxmissing=0.5, allow_stop_codons = False):
	# prune alignment 
	aln.subsetSamples(samples)
	# remove trailing stop codons
	aln.removeTrailingStopCodon()
	# filter codons with too many gaps
	try:
		if allow_stop_codons:
			break_on_stops = False
		else:
			break_on_stops = True
		aln.codonAwareFiltration(max_missingness = maxmissing, break_on_stops = break_on_stops)
	except:
		return None
	return aln

def run_codeml(alignment,tree,outfile,workdir, options = {},
	noisy = 3, verbose = 1, seqtype = 1, ndata = 1, icode =0,
	cleandata = 0, model = 0, NSsites = [0], CodonFreq = 7,
	clock = 0, fix_omega = 0, omega = 0.5
	):
	# initiate codeml object
	cml = codeml.Codeml(alignment = alignment, tree = tree, out_file = outfile, working_dir = workdir)
	# loop through the function arguments and set them
	cml.set_options(
		noisy = noisy,
		verbose = verbose,
		seqtype = seqtype,
		ndata = ndata,
		icode = icode,
		cleandata = cleandata,
		model = model,
		NSsites = NSsites,
		CodonFreq = CodonFreq,
		clock = clock,
		fix_omega = fix_omega,
		omega = omega
		)
	for arg in list(locals().keys())[5:-1]:
		cml.set_options(**{arg: eval(arg)})
	# then set any custom options
	if len(options.keys()) > 0:
		for key in options.keys():
			cml.set_options(**{key:options[key]})
	#		cml.set_options(key = options[key])
	#	cml.set_options(**arg = eval(arg))
	# run it
	results = cml.run(verbose=True)
	# return the codeml object
	#cml.print_options()
	return results

# define an argument parser
parser = argparse.ArgumentParser(description='Run branch site model in codeml to a null model using a likelihood ratio test')
# add arguments
parser.add_argument('-a', '--alignment', help='alignment file', required=True)
parser.add_argument('-t', '--tree', help='tree file', required=True)
parser.add_argument('-s', '--samples', help='samples to keep in analyses', required=False)
parser.add_argument('-w', '--workdir', help='working directory', default='./')
parser.add_argument('-n', '--name', help='name/prefix of analysis', default="codeml")
parser.add_argument('-f','--foreground', help='foreground clades', nargs='+', required=True)
parser.add_argument('-o', '--outgroup', help='outgroup', required=False, default=None)
# missingness filter
parser.add_argument('-m', '--maxmissing', help='maximum missingness allowed', default=0.5)
# degrees of freedom
parser.add_argument('-d', '--df', help='degrees of freedom', default=1)
# and all the codeml options
parser.add_argument('-noisy', help='verbose output', default=3)
parser.add_argument('-verbose', help='verbose output', default=1)
parser.add_argument('-seqtype', help='verbose output', default=1)
parser.add_argument('-ndata', help='verbose output', default=1)
parser.add_argument('-icode', help='verbose output', default=0)
parser.add_argument('-cleandata', help='verbose output', default=0)
parser.add_argument('-model', help='verbose output', default=0)
parser.add_argument('-NSsites', help='verbose output', default=[0])
parser.add_argument('-CodonFreq', help='verbose output', default=7)
parser.add_argument('-clock', help='verbose output', default=0)
parser.add_argument('-fix_omega', help='verbose output', default=0)
parser.add_argument('-omega', help='verbose output', default=0.5)
parser.add_argument('--allow-stop-codons', help='allow stop codons in alignment', action='store_true', default=False)

test = False # manually set to True if testing

# parse the arguments
args = parser.parse_args()

# if test == True:
# 	args = testargs(aln, tree, samples, outgroup, foreground, maxmissing, allow_stop_codons, workdir, name, df)

# start with reading the tree
tree = Tree(args.tree)

# if there are no sample argument, fetch the leaf names from the tree
if args.samples is None:
	samples = [i.name for i in tree.get_leaves()]
else:
	# otherwise read the samples file
	samples = []
	with open(args.samples, "r") as f:
		for line in f:
			samples.append(line.strip())

# make a list of the foreground branches from args
foreground = []
for clade in args.foreground:
	foreground.append(clade.split(","))



# prune and prep the tree
tree.prune(samples)
treestring = prep_tree(tree, foreground, args.outgroup)
# write the tree to a file in the working directory
with open(os.path.join(args.workdir, "{}_codeml_tree.nwk".format(args.name)), "w") as f:
	f.write(treestring)

# read the alignment
aln = fn.readSequenceFile(args.alignment)
# prep the alignment
aln = prep_alignment(aln, samples, float(args.maxmissing), allow_stop_codons=args.allow_stop_codons)
# if the alignment is None, exit
if aln is None:
	sys.exit("Check alignment: {}, probably internal stop codons or something messing it up.".format(args.alignment))
# and write it to a phylip file
fn.writeSequenceFile(aln, os.path.join(args.workdir, "{}_codeml_alignment.phy".format(args.name)))

# get the full paths to the alignment and tree files
alignment = os.path.abspath(os.path.join(args.workdir, "{}_codeml_alignment.phy".format(args.name)))
tree = os.path.abspath(os.path.join(args.workdir, "{}_codeml_tree.nwk".format(args.name)))

# specify options for null model
nulloptions = {
	'fix_omega':1, 
	'omega': 1, 
	'NSsites': [2], 
	'model':2
	}
# store the current path
cwd = os.getcwd()

# make a directory for the null model
os.mkdir(os.path.join(args.workdir, "nullmodel"))
# and change to it
os.chdir(os.path.join(args.workdir, "nullmodel"))

# run codeml null model
nullmodel = run_codeml(alignment = alignment, tree = tree, verbose = 1, outfile = "{}_codeml_null.out".format(args.name), workdir = "./", options = nulloptions)

# go back to initial working directory
os.chdir(cwd)

# make a directory for the alternative model
os.mkdir(os.path.join(args.workdir, "altmodel"))
# and change to it
os.chdir(os.path.join(args.workdir, "altmodel"))

# options for alternative model
altoptions = {
	'model':2, 
	'NSsites': [2], 
	'omega': 0.5,
	'fix_omega':0,
	}

# run codeml alternative model
alternative = run_codeml(alignment = alignment, tree = tree, outfile = "{}_codeml_alt.out".format(args.name), workdir = "./", options = altoptions)

# when these runs are done, we want to change back to the cwd variable
os.chdir(cwd)


# get the likelihoods
null_lnL = nullmodel.get("NSsites").get(2).get("lnL")
alt_lnL = alternative.get("NSsites").get(2).get("lnL")

# specify degrees of freedom
df = args.df

# calculate the likelihood ratio
LRT = 2*(alt_lnL - null_lnL)
# and the p-value if LRT is positive
if LRT > 0:
	p = cdf_chi2(df,LRT)
else:
	p = "NA"

# write the results to a file
with open(os.path.join(args.workdir, "{}_codeml_lrt-results.txt".format(args.name)), "w") as f:
	f.write("null_lnL\t{}\n".format(null_lnL))
	f.write("alt_lnL\t{}\n".format(alt_lnL))
	f.write("LRT\t{}\n".format(LRT))
	f.write("p\t{}\n".format(p))

# if the p-value is significant, we want to get the foreground branches
#if p < 0.05:
#	sys.stdout.write("p-value is significant! getting foreground branches\n")
#	# get the foreground branches
#	foreground_branches = alternative.get("NSsites").get(2).get("parameters").get("foreground w").keys()
#	# and write them to a file
#	with open(os.path.join(args.workdir, "{}_codeml_foreground_branches.txt".format(run_name)), "w") as f:
#		for branch in foreground_branches:
#			f.write("{}\n".format(branch))