#!/bin/bash

# This workflow will take the supercontig output of HybPiper and return a supercontig that
# contains heterozygous positions as ambiguity bases. Uses paired reads.

#The script should be run on a FASTA file containing all the supercontigs of interest.

if [[ $# -eq 0 ]] ; then
    echo 'usage: map_to_supercontigs.sh prefix readfile1.fq readfile2.fq'
    exit 1
fi

#############COMMAND LINE ARGUMENTS############

prefix=$1
read1fq=$2
read2fq=$3

supercontig=$prefix.supercontigs.fasta

#####STEP ZERO: Make Reference Databases

picard CreateSequenceDictionary \
R=$supercontig
bwa index $supercontig
samtools faidx $supercontig

#####STEP ONE: Map reads

echo "Mapping Reads"

bwa mem $supercontig $read1fq $read2fq | samtools view -bS - | samtools sort - -o $supercontig.sorted.bam

picard FastqToSam  \
F1=$read1fq \
F2=$read2fq \
O=$supercontig.unmapped.bam \
SM=$supercontig

picard MergeBamAlignment \
ALIGNED=$supercontig.sorted.bam \
UNMAPPED=$supercontig.unmapped.bam \
O=$supercontig.merged.bam \
R=$supercontig

#####STEP TWO: Mark duplicates

echo "Marking Duplicates"
picard MarkDuplicates \
I=$supercontig.merged.bam \
O=$supercontig.marked.bam \
M=$supercontig.metrics.txt

#######STEP THREE: Identify variants, select only SNPs

echo "Identifying variants"

samtools index $supercontig.marked.bam

gatk HaplotypeCaller \
-R $supercontig \
-I $supercontig.marked.bam \
-O $supercontig.vcf

time gatk SelectVariants \
-R $supercontig \
-V $supercontig.vcf \
-select-type SNP \
-O $supercontig.snps.vcf

######STEP FOUR: Output new supercontig FASTA with ambiguity codes

echo "Generating IUPAC FASTA file"

gatk FastaAlternateReferenceMaker \
-R $supercontig \
-O $supercontig.iupac.fasta \
-V $supercontig.snps.vcf \
--use-iupac-sample $supercontig
