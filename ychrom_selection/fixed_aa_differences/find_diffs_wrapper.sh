# this is a wrapper around the fixedAaDiffs.py script

#### Set parameters

# popfile to use
popfile=mitis-denti_vs_pog-wolfi.popfile.txt

# directory with the alignments
FASTADIR=../gene_alignments/

OUTDIR=AA_FIXED_DIFFS
mkdir -p $OUTDIR
mkdir -p ${OUTDIR}/aa_alignments

NAME=mitden_vs_pogwol

#####

echo "" > ${OUTDIR}/${NAME}.log

# loop through the alignments
for file in $FASTADIR/*.phy
do
# fetch the gene name
gene=$(basename ${file} | cut -d "_" -f 1)
echo $gene >> ${OUTDIR}/${NAME}.log

# run the fixedAaDiffs.py script
python fixedAaDiffs.py -a $file -p $popfile -o ${OUTDIR}/${NAME}.${gene} --write-aln ${OUTDIR}/aa_alignments/${NAME}.${gene}_aa.fa \
--translate >> ${OUTDIR}/${NAME}.log
done

# parse counts
echo -e "gene\tgoodsites\tfixed_diffs\tmissing_sites" > ${OUTDIR}/${NAME}.all_counts.txt

# make a table with the number of fixed differences per gene
for file in ${OUTDIR}/${NAME}*fixed_diffs_counts.txt
do
gene=$(basename ${file} | cut -d "." -f 2 | cut -d "_" -f 1)
goodsites=$(tail -n +2 $file | cut -f 1)
fixed=$(tail -n +2 $file | cut -f 2)
missing=$(tail -n +2 $file | cut -f 3)
echo -e "${gene}\t${goodsites}\t${fixed}\t${missing}" >> ${OUTDIR}/${NAME}.all_counts.txt
done

# some python lines to check samples that have internal stops for each gene
python3 parse_logfile.py ${OUTDIR}/${NAME}
