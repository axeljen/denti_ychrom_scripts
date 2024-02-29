# this script will combine the outputs from chromosome-wise runs of dsuite dtrios. Just give the dir with all the output, the tree used in the run and a name for the output as command line inputs.

# directory to combine as first input
INDIR=$1
INDIR=$(realpath ${INDIR})
# tree to use when combining as second
TREE=$2
TREE=$(realpath ${TREE})
# name to give to the output
NAME=$3

# make a list of the combine files 
files=$(find ${INDIR} -name "*_Dmin.txt" | grep -v -E "NC_023671|NC_041774" | rev | cut -d "_" -f 2- | rev)
cd ${INDIR}
Dsuite DtriosCombine -o ${NAME} -t ${TREE} ${files}