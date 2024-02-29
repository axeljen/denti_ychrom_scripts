#!/bin/bash

#SBATCH -A snic2022-5-561
#SBATCH -p node
#SBATCH -N 1
#SBATCH -C mem512GB
#SBATCH -t 2-00:00:00
#SBATCH -J mitofinder
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=axel.jensen@ebc.uu.se
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

#we'll use bbmap for subsampling input reads, and trimmomatic to trim those, load these modules
module load bioinfo-tools bbmap trimmomatic samtools

#mitofinder is installed manually, specify the path to your installation here NOTE: you need to change this to your own installation
MITOFINDER_PATH=/domus/h1/axeljen/MitoFinder

#specify number of threads
THREADS=16

ASSEMBLER="--metaspades"

#specify the path to the reference genome to work on, instructions on how to find and download this is available on mitofinder's github
REFERENCE=/crex/proj/sllstore2017021/nobackup/GUENON_GENOMES/REFERENCE/MITOGENOME/MITOFINDER_REFERENCE/chl.sab_NC_008066.1.mito.gb

# take a bam/cram file as input
BAM=$1

OUTDIR=$2

# grab sample name from this
SAMPLE=$(samtools view -H ${BAM} | grep @RG | cut -f 3 | cut -d ":" -f 2 | head -n 1)

mkdir -p ${OUTDIR}/${SAMPLE}

OUT=$(realpath ${OUTDIR}/${SAMPLE})

# to extract paired reads from the cramfile, first need to sort it by read
samtools sort -n -@ ${THREADS} -o ${SNIC_TMP}/sort.${SAMPLE}.bam ${BAM}

# then extract paired fq reads
samtools fastq -1 ${SNIC_TMP}/${SAMPLE}_1.fastq -2 ${SNIC_TMP}/${SAMPLE}_2.fastq \
	-0 /dev/null -s /dev/null ${SNIC_TMP}/sort.${SAMPLE}.bam 

# set these variables for easier parsing
READS_1=${SNIC_TMP}/${SAMPLE}_1.fastq
READS_2=${SNIC_TMP}/${SAMPLE}_2.fastq

#first thing we'll subsample the reads as mitochondrial reads are usually overrepresented in short read data, and using all data is unnecessary slow and may actually give worse results.
NREADS=3000000 # can be worth bumping this if the assembly is poor

#subsampling using bbmap, putting the subsampled reads on ${SNIC_TMP}
reformat.sh t=${THREADS} \
	in1=${READS_1} in2=${READS_2} \
	out1=${SNIC_TMP}/sub.$(basename ${READS_1}) out2=${SNIC_TMP}/sub.$(basename ${READS_2}) \
	samplereadstarget=${NREADS}

#now we trim adapters off of the subsampled reads, still working on tempstorage
trimmomatic PE -threads ${THREADS} \
 	${SNIC_TMP}/sub.$(basename ${READS_1}) ${SNIC_TMP}/sub.$(basename ${READS_2}) \
	${SNIC_TMP}/tr.paired.$(basename ${READS_1}) ${SNIC_TMP}/tr.u.$(basename ${READS_1}) \
	${SNIC_TMP}/tr.paired.$(basename ${READS_2}) ${SNIC_TMP}/tr.u.$(basename ${READS_2}) \
	ILLUMINACLIP:${TRIMMOMATIC_ROOT}/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36


#let's assemble on tempstorage, and then fetch the output when finished
cd ${SNIC_TMP}

#now let's run mitofinder, the only option I added from the default is the --adjust-direction, which makes sure that all mitochondrias are outputted in the same direction
${MITOFINDER_PATH}/mitofinder --adjust-direction \
	${ASSEMBLER} \
	-t mitfi \
	-r ${REFERENCE} \
	-j ${SAMPLE} \
	-1 ${SNIC_TMP}/tr.paired.$(basename ${READS_1}) \
	-2 ${SNIC_TMP}/tr.paired.$(basename ${READS_2}) \
	-o 2 -p ${THREADS}

# fetch all the output from the Final_Results dir, and leave the rest behind
cp ${SAMPLE}/${SAMPLE}*MitoFinder*Final_Results*/* ${OUT}

# get some info from this sample
# number of countigs
contigs_found=$(find ${OUT} -name "*_contig*.fasta" | grep -v -E "_NT|_AA" | wc -l)

# coverages, lengths and circularization for these
for i in $(seq 1 ${contigs_found})
do
# if single contig, we just grab this ones coverage
if [[ ${contigs_found} == 1 ]];
then
name=$(basename $(find ${OUT} -name "*_contig.fasta"))
cov=$(cat ${OUT}/*.infos | sed -n 1p | rev | cut -d "_" -f 1 | rev)
length=$(cat ${OUT}/*.infos | sed -n 1p | rev | cut -d "_" -f 3 | rev)
echo "${SAMPLE};${name};${cov};${length};$(basename ${REFERENCE});$(date | sed "s| |_|g" | sed "s|:|-|g");${NREADS}" >> $(dirname ${OUT})/assembly_stats.txt
else
name=$(basename $(find ${OUT} -name "*_contig_${i}.fasta"))
cov=$(cat ${OUT}/*contig_${i}.infos | sed -n 3p | rev | cut -d "_" -f 1 | rev)
length=$(cat ${OUT}/*contig_${i}.infos | sed -n 3p | rev | cut -d "_" -f 3 | rev)
echo "${SAMPLE};${name};${cov};${length};$(basename ${REFERENCE});$(date | sed "s| |_|g" | sed "s|:|-|g");${NREADS}" >> $(dirname ${OUT})/assembly_stats.txt
fi
done

