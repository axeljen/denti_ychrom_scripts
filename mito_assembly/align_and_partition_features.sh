# wrapper script to concatenate the extracted features/codonbases into different blocks and then a large alignment with a partitionfile accompanying. 
# Used as input for the iqtree partitioned phylogeny.

# prepared a list with all features to concatenate here:
FEATURES=featurelist.txt

# loop through all of them and find all corresponding phylips, first double checking so that everything adds up numberwise
for f in $(cat ${FEATURES})
do
echo ${f}
find extracted_features -name "*_${f}_*" | grep -v full_genes | wc -l
done
module load bioinfo-tools MAFFT

# start with concatenating gene codons, as these needs some extra care 
for f in $(cat ${FEATURES} | head -n 13)
do
echo ${f}
for i in first second third
do
for phylip in $(find extracted_features/genes/${i}_codons -name "*_${f}_*")
do 
name=$(cat ${phylip} | tail -n +2 | cut -d " " -f 1)
seq=$(cat ${phylip} | tail -n +2 | awk ' { print $2 } ')
echo ">${name}" >> ${f}_${i}_codon_sequences.fa
echo ${seq} >> ${f}_${i}_codon_sequences.fa
done
# align with mafft
mafft ${f}_${i}_codon_sequences.fa > ${f}_${i}_codon_aligned.fa
done
done

# moved them all to be in their respecive subdir to extracted_features/genes/*codon_position

# next, we'll concatenate all the trnas for all the samples
for i in $(find extracted_features/tRNA/trnas_concat -name "*.phy")
do
name=$(cat ${i} | tail -n +2 | cut -d " " -f 1)
seq=$(cat ${i} | tail -n +2 | awk ' { print $2 } ')
echo ">${name}" >> extracted_features/tRNA/trnas_concat/trnas_concat_sequences.fa
echo ${seq} >> extracted_features/tRNA/trnas_concat/trnas_concat_sequences.fa
done
# align these too with mafft
mafft extracted_features/tRNA/trnas_concat/trnas_concat_sequences.fa > extracted_features/tRNA/trnas_concat/trnas_concat_aligned.fa

# and last the rRNA
for f in rrnL rrnS
do
for i in $(find extracted_features/rRNA/${f} -name "*.phy")
do
name=$(cat ${i} | tail -n +2 | cut -d " " -f 1)
seq=$(cat ${i} | tail -n +2 | awk ' { print $2 } ')
echo ">${name}" >> extracted_features/rRNA/${f}/${f}_concat_sequences.fa
echo ${seq} >> extracted_features/rRNA/${f}/${f}_concat_sequences.fa
done
mafft extracted_features/rRNA/${f}/${f}_concat_sequences.fa > extracted_features/rRNA/${f}/${f}_concat_aligned.fa
done

# now, make a list with all the aligned files that we will concatenate into one
find extracted_features -name "*aligned.fa" > alignments_to_concat.txt

# and then use my python script for concatenating
python3 ~/phylogenomics/concatSequences.py -i alignments_to_concat.txt -o trimmed_concatenated_features_aligned.phy --raxml-partitions