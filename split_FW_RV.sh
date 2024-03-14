#!/bin/bash

# $1 = .._1.fastq
# $2 = .._2.fastq
# $3 = primer_file

touch log.txt
touch report.txt

# INPUT FILES

zcat -f "$1" > input1.fastq
echo -n "input1 " >> report.txt && wc -l input1.fastq | awk '{x=$1/4; print x}' >> report.txt

zcat -f "$2" > input2.fastq
echo -n "input2 " >> report.txt && wc -l input2.fastq | awk '{x=$1/4; print x}' >> report.txt


# ADAPTOR TRIMMING

cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -m 60 -q 18 -j 4 -o file_1_trimmed.fastq -p file_2_trimmed.fastq input1.fastq input2.fastq >> log.txt
echo -n "file_1_trimmed " >> report.txt && wc -l file_1_trimmed.fastq | awk '{x=$1/4; print x}' >> report.txt
echo -n "file_2_trimmed " >> report.txt && wc -l file_2_trimmed.fastq | awk '{x=$1/4; print x}' >> report.txt

rm input1.fastq
rm input2.fastq

# SEPARATE FORWARD AND REVERSE READS

# Split the forward and reverse reads in each of the two files 
cat file_1_trimmed.fastq | fastx_barcode_splitter.pl --bcfile "$3" --bol --prefix file_1_trimmed- --suffix .fastq --mismatches 5 --partial 5
cat file_2_trimmed.fastq | fastx_barcode_splitter.pl --bcfile "$3" --bol --prefix file_2_trimmed- --suffix .fastq --mismatches 5 --partial 5

# Concatenate FORWARD reads from file1 and file2. Contatenate REVERSE reads from file1 and file2.
cat file_1_trimmed-FW.fastq file_2_trimmed-FW.fastq > FWreads.fastq
cat file_1_trimmed-RV.fastq file_2_trimmed-RV.fastq > RVreads.fastq

# Remove the partial FW and RV only files
rm file_1_trimmed-FW.fastq
rm file_1_trimmed-RV.fastq
rm file_2_trimmed-FW.fastq
rm file_2_trimmed-RV.fastq

fastqc FWreads.fastq
fastqc RVreads.fastq


