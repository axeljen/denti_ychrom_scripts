# wrapper script to fetch the annotated features from the mitogenomes to separate sequences, and split the pcgs by codon position

# assumes that the annotated mitochondria are in a dir called annotated_assemblies
# and that there's an output dir called extracted_features

# will loop through all fasta sequences in annotated_assemblies dir
for fasta in $(find annotated_assemblies -name "*.fa" -mindepth 1 -maxdepth 1)
do
# fetch the corresponding gff file
gff=${fasta/.fa/.gff}
# use the subsetSequence script to fetch CDS, tRNA and rRNA
python3 /domus/h1/axeljen/phylogenomics/subsetSequence.py -i ${fasta} -o extracted_features/$(basename ${fasta/.fa}) \
	-gff ${gff} \
	--extract-type CDS,tRNA,rRNA

# next, move the alignments to their respective parent directories

## starting with the tRNAs
mv extracted_features/$(basename ${fasta/.fa})_tRNA-*.phy extracted_features/tRNA/all_trnas
# concatenate all tRNAs into a single sequence, will just use the first gff file to set the order, to make sure they're all concatenated in the same way
cat $(find annotated_assemblies -name "*.gff" | head -n 1) | awk ' $3 == "tRNA" { print $9 } ' | cut -d "=" -f 2 > trna_concat_order.txt
# then find the correct tRNA's for this sample and concat
echo "" > $(basename ${fasta/.fa/}).trnas_concat.txt
for trna in $(cat trna_concat_order.txt)
do
find extracted_features/tRNA/all_trnas -name "$(basename ${fasta/.fa/})*${trna}_*.phy" >> $(basename ${fasta/.fa/}).trnas_concat.txt
done
# concatenate them
python3 /domus/h1/axeljen/phylogenomics/concatSequences.py -i $(basename ${fasta/.fa/}).trnas_concat.txt -o extracted_features/tRNA/trnas_concat/$(basename ${fasta/.fa/})_tRNAs.phy
# and then work the rRNAs
mv extracted_features/$(basename ${fasta/.fa})_rrnS_*.phy extracted_features/rRNA/rrnS/
mv extracted_features/$(basename ${fasta/.fa})_rrnL_*.phy extracted_features/rRNA/rrnL/
# and last the genes, which should be the only ones remaining now
mv extracted_features/$(basename ${fasta/.fa})*_*.phy extracted_features/genes/full_genes
# go through all of them and chop into codons
for gene in extracted_features/genes/full_genes/$(basename ${fasta/.fa/}*.phy)
do
python3 /domus/h1/axeljen/phylogenomics/codonChopper.py -i ${gene} -o extracted_features/genes --strict-triplets --outformat phy
done
# move codons to their respective dirs
mv extracted_features/genes/*codon_1.phy extracted_features/genes/first_codons
mv extracted_features/genes/*codon_2.phy extracted_features/genes/second_codons
mv extracted_features/genes/*codon_3.phy extracted_features/genes/third_codons
done