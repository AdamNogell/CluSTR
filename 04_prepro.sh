#!/bin/bash

# $1 = reads1.fq.gz
# $2 = reads2.fq.gz
# $3 = primers.txt
# $4 = reference_genome.fasta

# CREATE NEEDED DIRECTORIES AND FILES
mkdir split
mkdir top
mkdir mapping
touch report.txt
touch log.txt
echo "Number of reads" >> report.txt

# UNZIP READS
zcat -f "$1" > input1.fastq
echo -n "input file 1: " >> report.txt &&
    wc -l input1.fastq | awk '{print int($1/4)}' >> report.txt

zcat -f "$2" > input2.fastq
echo -n "input file 2: " >> report.txt &&
    wc -l input2.fastq | awk '{print int($1/4)}' >> report.txt

cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -m 60 -q 18 -j 4 -o split/file_1_trimmed.fastq -p split/file_2_trimmed.fastq input1.fastq input2.fastq >> log.txt
echo -n "file 1 after adapter cutting: " >> report.txt && 
    wc -l split/file_1_trimmed.fastq | awk '{print int($1/4)}' >> report.txt
echo -n "file 2 after adapter cutting: " >> report.txt && 
    wc -l split/file_2_trimmed.fastq | awk '{print int($1/4)}' >> report.txt

fastx_barcode_splitter.pl --bcfile "$3" --bol --prefix split/file_1_trimmed- --suffix .fastq --mismatches 5 --partial 5 < split/file_1_trimmed.fastq >> log.txt
echo -n "file 1 FW reads: " >> report.txt && 
    wc -l split/file_1_trimmed-FW.fastq | awk '{print int($1/4)}' >> report.txt
echo -n "file 1 RV reads: " >> report.txt && 
    wc -l split/file_1_trimmed-RV.fastq | awk '{print int($1/4)}' >> report.txt
fastx_barcode_splitter.pl --bcfile "$3" --bol --prefix split/file_2_trimmed- --suffix .fastq --mismatches 5 --partial 5 < split/file_2_trimmed.fastq >> log.txt
echo -n "file 2 FW reads: " >> report.txt && 
    wc -l split/file_2_trimmed-FW.fastq | awk '{print int($1/4)}' >> report.txt
echo -n "file 2 RV reads: " >> report.txt && 
    wc -l split/file_2_trimmed-RV.fastq | awk '{print int($1/4)}' >> report.txt

cat split/file_1_trimmed-FW.fastq split/file_2_trimmed-FW.fastq > split/FWreads.fastq
echo -n "FW reads total: " >> report.txt && 
    wc -l split/FWreads.fastq | awk '{print int($1/4)}' >> report.txt
cat split/file_1_trimmed-RV.fastq split/file_2_trimmed-RV.fastq > split/RVreads.fastq
echo -n "RV reads total: " >> report.txt && 
    wc -l split/RVreads.fastq | awk '{print int($1/4)}' >> report.txt

rm input1.fastq
rm input2.fastq
rm split/file_1_trimmed-FW.fastq
rm split/file_1_trimmed-RV.fastq
rm split/file_2_trimmed-FW.fastq
rm split/file_2_trimmed-RV.fastq

# FASTQC - UNHASH IF NEEDED
#fastqc FWreads.fastq
#fastqc RVreads.fastq

cutadapt -g GCAGAATCGGTAGGCTCTTC -j 4 -o top/RVreads_primerless.fastq split/RVreads.fastq >> log.txt
echo -n "RV reads after primer cut: " >> report.txt && 
    wc -l top/RVreads_primerless.fastq | awk '{print int($1/4)}' >> report.txt
awk '{NR%4 == 2}' top/RVreads_primerless.fastq > top/RVreads_primerless_awk.txt
sort top/RVreads_primerless_awk.txt | uniq -c | awk '{print $1 "\t" $2}' | sort -nr -k 1 | head -n 200 > top/RVreads_top200.tsv

# MAPPING
bwa index "$4" >> log.txt
bwa mem -t 4 "$4" top/RVreads_primerless.fastq > mapping/RVreads.sam
samtools view -b mapping/RVreads.sam > mapping/RVreads.bam
samtools sort mapping/RVreads.bam -o mapping/RVreads_sorted.bam
samtools index mapping/RVreads_sorted.bam
samtools view -b -o mapping/RV_region_sorted.bam mapping/RVreads_sorted.bam U13369.1:41469-41789
samtools index mapping/RV_region_sorted.bam