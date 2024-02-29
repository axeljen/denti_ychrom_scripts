# this is a wrapper around the fixedAaDiffs.py script

# popfile to use
popfile=mitis-denti_vs_pog-wolfi.popfile.txt

# directory with the fasta files
FASTADIR=Y_TRANSCRIPT_ALIGNMENTS

NAME=mitden_vs_pogwol_restr
echo "" > AA_FIXED_DIFFS/${NAME}.log

# loop through the fasta files
for file in $FASTADIR/*.fa
do
# fetch the gene name
gene=$(basename ${file} | cut -d "_" -f 1)
echo $gene >> AA_FIXED_DIFFS/${NAME}.log
# first a small python script to make sure that both denti AND at least one other mitis is not missing more than 10% of the sites
# run the script
python fixedAaDiffs.py -a $file -p $popfile -o AA_FIXED_DIFFS/${NAME}.${gene} --write-aln AA_ALIGNMENTS/${NAME}.${gene}_aa.fa \
--translate >> AA_FIXED_DIFFS/${NAME}.log
done

echo -e "gene\tgoodsites\tfixed_diffs\tmissing_sites" > AA_FIXED_DIFFS/${NAME}.all_counts.txt

# make a table with the number of fixed differences per gene
for file in AA_FIXED_DIFFS/${NAME}*fixed_diffs_counts.txt
do
gene=$(basename ${file} | cut -d "." -f 2 | cut -d "_" -f 1)
goodsites=$(tail -n +2 $file | cut -f 1)
fixed=$(tail -n +2 $file | cut -f 2)
missing=$(tail -n +2 $file | cut -f 3)
echo -e "${gene}\t${goodsites}\t${fixed}\t${missing}" >> AA_FIXED_DIFFS/${NAME}.all_counts.txt
done

# some python lines to check samples that have internal stops for each gene
python3 parse_logfile.py AA_FIXED_DIFFS/${NAME}
