
# script assumes we're in a directory containing all the mitofinder assemblies from previous step (mitofinder_asm.sh)
mkdir RAW # we'll store the raw mitofinder assemblies here
mv *.fasta RAW
mkdir -p ROT # and the rotated ones here

# loop through all the mitofinder assemblies and rotate them to the same startpos
for i in $(find RAW -name *.fasta)
do
python3 ~/phylogenomics/circularizeAndRotate.py -i ${i} -o ROT/rot.$(basename ${i}) --force-rotation -ref /crex/proj/sllstore2017021/nobackup/DENTI_LESULA_PROJECT/ANALYSES/MITOGENOMES_NEW/chlsab_ref_NC_008066.1.fa
done

# now it should be fairly straightforward to align those and then trim manually
module load bioinfo-tools MAFFT

cd ROT
cat *.fasta > mitogenomes_raw_rotated.fa

# run mafft on this
mafft mitogenomes_raw_rotated.fa > aln.mitogenomes_raw_rotated.fa

# I then used this alignment to manually trim any erroneous ends/repeats in the mitogenomes, re-extract each sequence and then we run mitofinder_annotate.sh on each of them