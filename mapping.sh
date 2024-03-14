#!/bin/bash

# $1 = reference_genome.fasta
# $2 = FWreads.fastq
# $3 = region_of_interest.txt

bwa index "$1"

bwa mem -t 4 "$1" "$2" > RVreads.sam

samtools view -b RVreads.sam > RVreads.bam

samtools sort RVreads.bam -o RVreads_sorted.bam

samtools index RVreads_sorted.bam

samtools view -b -o region_sorted_output.bam RVreads_sorted.bam "$3"

samtools index region_sorted_output.bam
