# This file contains the steps to sample random regions from a reference genome, and create alignments for each region with the prespecified samples.
# The last steps generates an alignment file in the bpp input format.
# Note that these steps should be run interactively, one step at a time, since some of them need manual evaluation/intervention.

## File paths:
# path to directory with filtered vcf files. I had one vcf file per chromosome, so this may need modification if you have a different setup.
VCFDIR=/path/to/vcf/files
# path to the reference genome
REFERENCE=/path/to/reference/genome
# path to a file with the names of the samples to include in the alignment
SAMPLES=samples_in_bpp.txt
# path to a file with the regions to sample from (I made a file with regions located at least 10kb from the nearest gene in the macaque reference genome)
REGIONS=noncoding_regions_autosomes.txt

## Number of loci to collect
NLOCI=1000

# path to store the alignments
OUTDIR=noncoding_loci
mkdir -p ${OUTDIR}
for i in $(seq 1 ${NLOCI})
do
echo "Starting with locus ${i}..."
while :
do
# fetch a random interval/region
randint=$(shuf -i 1-$(cat ${REGIONS} | wc -l) -n 1)
region=$(cat ${REGIONS} | sed -n ${randint}p)
# fetch a random 1k window from this region
chrom=$(echo ${region} | cut -d " " -f 1)
min=$(echo ${region} | cut -d " " -f 2)
max=$(($(echo ${region} | cut -d " " -f 3) - 1000))
start=$(shuf -i ${min}-${max} -n 1)
end=$(( ${start} + 1000 -1))
# get the vcf file for this chromosome, again this may need modification with a different setup
VCF=$(find ${VCFDIR} -name "*${chrom}*.vcf.gz")
# make an alignment
python3 ~/phylogenomics/vcfToMSA.py -v ${VCF} -o ${OUTDIR}/${chrom}_${start}_${end}.phy -of phylip --haploidize IUPAC -r ${chrom}:${start}-${end} --samples ${SAMPLES} \
--reference ${REFERENCE}
# get the sequence length after removing all the NS
filtered_seqlen=$(python3 check_filtered_len.py ${OUTDIR}/${chrom}_${start}_${end}.phy)
if [ $filtered_seqlen -gt 799 ]
then
# if the alignment is long enough, break the loop and move on to the next locus
echo "Finished with locus ${i}"
break
else
# Discarding loci where more than 200 bp have N's in them, and trying again
echo "This alignment discarded, trying out a new one..."
rm -f ${OUTDIR}/${chrom}_${start}_${end}.phy
fi
done
done
done

# this step checks that loci are at least 50 k apart from each other, to avoid linkage.
last_end=0
for i in $(ls -1 ${OUTDIR} | sort -V)
do
start=$(echo ${i} | cut -d "_" -f 3)
end=$(echo ${i} | cut -d "_" -f 4 | cut -d "." -f 1)
echo ${i} $(( ${start} - ${last_end} ))
last_end=${end} 
done | awk ' $2 < 50000 && $2 > 0 ' | cut -f 1 -d " " > to_remove.txt


# remove those and resample until we have 1000 that are at least 50 k apart ( set NLOCI above to the number of lines in to_remove.txt, and rerun the loop )
# this needs to be repeated until the previous step returns an empty file (to_remove.txt) 
for phy in $(cat to_remove.txt)
do
rm -f noncoding_1/${phy}
done

# remove the to_remove.txt file
rm to_remove.txt

# convert to bpp format, one file with all loci where all samples are prefixed by ^
> ${OUTDIR}/1000_loci_bpp_format.txt
# manually specify number of samples and length of loci, since this will be the same for all loci
# In my case I had 9 samples and 1000 bp long loci
NSAMPLES=9
LOCILENGTH=1000
for i in $(ls -1 ${OUTDIR} | sort -V | grep -v loci_bpp_format.txt)
do
echo -e "\t${NSAMPLES}\t${LOCILENGTH}">> ${OUTDIR}/1000_loci_bpp_format.txt
while read line
do
echo "^${line}" >> ${OUTDIR}/1000_loci_bpp_format.txt
done < <(tail -n +2 ${OUTDIR}/${i})
done

# the vcfToMSA.py script adds the reference sequence under the name "reference" to the alignment. Since I don't want that (I have the macaque already included as another sample),
# I remove the reference sequence from the alignment file
cat ${OUTDIR}/1000_loci_bpp_format.txt | grep -v reference > tmp && mv tmp ${OUTDIR}/1000_loci_bpp_format.txt

# clean up
mv ${OUTDIR}/1000_loci_bpp_format.txt 1000_loci_bpp_format.txt
rm -r ${OUTDIR}