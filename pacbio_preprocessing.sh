#!/bin/bash

# pacbio_preprocessing.sh
# Author: Kirsten Reimer
# Sep 2019
# This script processes raw FASTQ PacBio reads and maps processed reads to the genome.
# Usage: pacbio_preprocessing.sh raw_file_pc.fastq sampleID

################################################################################
# This fist step must be done on my computer, since porechop is not installed on
# the Farnam cluster. All other steps in the script should be submitted as a
# batch job on the cluster.

# porechop -i raw_file.fastq -o raw_file_pc.fastq --extra_end_trim 0 --extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0 --min_split_read_size 100

################################################################################
#load required modules on Farnam cluster
  module load cutadapt PRINSEQ minimap2 SAMtools BEDTools

  FASTQ=$1
  ID=$2
  GENOME=/gpfs/ysm/datasets/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa

# make a new folder for each raw file to keep intermediates and fully processed reads
  cd /home/kar86/project/pacbio_data
  mkdir /home/kar86/project/pacbio_data/"$ID"
  cp /home/kar86/project/pacbio_data/"$FASTQ" /home/kar86/project/pacbio_data/"$ID"
  cd /home/kar86/project/pacbio_data/"$ID"

# cutadapt filters for the 3'end DNA adapter and trims it if it is found
# -n 1 since it should be found only once per read since chimeras are already split up using Porechop
# allowed error rate of 10% in matching sequence and minimum length 15 required to keep a trimmed read
# keep untrimmed output as these are reads in the opposite orientation that need to be trimmed in a second step
  cutadapt -a CTGTAGGCACCATCAAT -e 0.1 -m 15 --untrimmed-output=untrimmed.fastq --quiet -o trim1.fastq $FASTQ

# make reverse complement of untrimmed reads so they can be filtered again
  cat untrimmed.fastq | tr '\n' '	' | sed 's/@m1/\
@m1/g' | sed '/^$/d' > untrimmed.tab

  cut -f 1 untrimmed.tab > ID
  cut -f 2 untrimmed.tab | rev | sed 's/A/1/g' | sed 's/C/2/g' | sed 's/G/3/g' | sed 's/T/4/g' | sed 's/1/T/g' | sed 's/2/G/g' | sed 's/3/C/g' | sed 's/4/A/g' > SEQ
  cut -f 3 untrimmed.tab > PLUS
  cut -f 4 untrimmed.tab | rev > QUAL

  paste ID SEQ PLUS QUAL | sed 's/	/\
/g' > untrimmed_rev.fastq

  rm ID SEQ PLUS QUAL untrimmed.tab

# repeat cutadapt filtering on RC half of reads to get rid of 3'end adapter and discard untrimmed reads (anything not ligated to 3' end adapter)
  cutadapt -a CTGTAGGCACCATCAAT -e 0.1 -m 15 --discard-untrimmed --quiet -o trim2.fastq untrimmed_rev.fastq
  cat trim1.fastq trim2.fastq > "$ID"_adapters_trimmed.fastq

  rm trim1.fastq trim2.fastq untrimmed_rev.fastq untrimmed.fastq

# prinseq to remove PCR copies - exact duplicates and reverse complement duplicates
# then remove 5 nucleotides from left and right (XXXXX added by RT and UMI of 3'adapter)

  prinseq-lite.pl -fastq "$ID"_adapters_trimmed.fastq -derep 1 -out_good "$ID"_dedup -out_bad "$ID"_PCRduplicates
  prinseq-lite.pl -fastq "$ID"_dedup.fastq -trim_left 6 -trim_right 5 -out_good "$ID"_dedup_trimmed

# map to genome using minimap2, output as .sam file
  minimap2 -ax splice -uf -C5 --secondary=no $GENOME "$ID"_dedup_trimmed.fastq > "$ID".sam

# convert sam to bam, sort, and index, and convert to BED12 for downstream analysis
	samtools view -S -b "$ID".sam > tmp.bam
	samtools sort tmp.bam -o "$ID"_sorted.bam
	samtools index "$ID"_sorted.bam
	bamToBed -bed12 -i "$ID"_sorted.bam > "$ID"_sorted_12.bed
	rm tmp.bam

# move sam, bam, and bed files to folders for organization
  cp "$ID".sam /home/kar86/project/pacbio_data/sam_files
  cp "$ID"_sorted.bam /home/kar86/project/pacbio_data/bam_files
  cp "$ID"_sorted.bam.bai /home/kar86/project/pacbio_data/bam_files
  cp "$ID"_sorted_12.bed /home/kar86/project/pacbio_data/bed12_files
