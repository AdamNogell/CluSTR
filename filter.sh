#!/bin/bash

# filter all reads that were clipped from the left side and pass through the whole Q1
samtools view -H RVreads_sorted.bam > header1.sam && samtools view RVreads_sorted.bam | awk -v start=41667 -v end=41686 '{if ($4 <= start && $4 + length($10) >= end && $6 ~ /^[0-9]+S/) {split($6, arr, "S"); if ($4 <= start && ($4 + length($10) - arr[1]) >= end) print}}' | cat header.sam - | samtools view -b > filtered1.bam && rm header1.sam

# filter all reads that were not clipped from the left side and pass through the whole Q1
samtools view -H RVreads_sorted.bam > header2.sam && samtools view RVreads_sorted.bam | awk -v start=41667 -v end=41686 '{if ($4 <= start && $4 + length($10) >= end && $6 !~ /^[0-9]+S/) print}' | cat header.sam - | samtools view -b > filtered2.bam && rm header2.sam

# merge all reads that pass through Q1 and indexe the merged bam file
samtools merge final.bam filtered1.bam filtered2.bam && samtools index final.bam

#remove unnecessary files
rm filtered*



