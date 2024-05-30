import argparse
import sys
import random
# change to correct path to phylogenomics if needed
sys.path.insert(1, '/home/axeljen/phylogenomics')
import functions as fn

# parse the alignment from standard in
msa = fn.readSequenceFile(sys.argv[1])

# count the number of positions where any sequence has an 'N'
missing_count = 0
for i in range(0,msa.length):
    for sample in msa.sequences.keys():
        if msa.sequences[sample].sequence[i] == "N":
            missing_count += 1
            break

sys.stdout.write(str(msa.length - missing_count))
