import sys

## wrapper script to take the output from autosomal denti simulations, read the combined Dstat output,
## switch P1 and P2 and sign of Dstat if needed and add the gene flow proportion and simulation name to the row

# Dstat file, simname and gene flow proportion as arguments
dstatfile = sys.argv[1]
simname = sys.argv[2]
gfprop = sys.argv[3]

# read the Dstat file
dstat = open(dstatfile).readlines()

# zip to a dict
dstatdict = dict(zip(dstat[0].strip().split(), dstat[1].strip().split()))

# switch P1 and P2 if needed (if denti is P1, we'd want this represented as negative D instead, with wolfi still as P1)
if dstatdict['P1'] == 'denti':
	dstatdict['P1'] = 'wolfi'
	dstatdict['P2'] = 'denti'
	dstatdict['Dstatistic'] = str(-float(dstatdict['Dstatistic']))
	dstatdict['f4-ratio'] = str(-float(dstatdict['f4-ratio']))

# add the gene flow proportion and simulation name to the row
dstatdict['gfprop'] = gfprop
dstatdict['simname'] = simname

# print tab separated values to standard out
sys.stdout.write('\t'.join([dstatdict['simname'], dstatdict['gfprop'], dstatdict['P1'], dstatdict['P2'], dstatdict['P3'], dstatdict['Z-score'], dstatdict['Dstatistic'], dstatdict['f4-ratio']]) + '\n')
	