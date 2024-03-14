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

# PAIR MERGING

fastq-join -x file_1_trimmed.fastq file_2_trimmed.fastq -o file_joined.fastq
echo -n "merged " >> report.txt && wc -l file_joined.fastqjoin | awk '{x=$1/4; print x}' >> report.txt

# UNIFYING READ ORIENTATION

cat file_joined.fastqjoin | fastx_barcode_splitter.pl --bcfile "$3" --bol --prefix file_join_split- --suffix .fastq --mismatches 5 --partial 5
seqkit seq -r -p file_join_split-RV.fastq > file_join_split_RVrc.fastq
cat file_join_split-FW.fastq file_join_split_RVrc.fastq > sense.fastq

echo -n "unified " >> report.txt && wc -l sense.fastq | awk '{x=$1/4; print x}' >> report.txt

# PRIMER CLIPPING

fwprimer=$(cat "$3" | cut -f 2  | sed -n '1p')
rvprimer=$(cat "$3" | cut -f 2  | sed -n '2p' | rev | tr ATGC TACG) # This primer now also needs to be a revcom

cutadapt -g "$fwprimer" -j 4 -o sense_FW_clipped.fastq sense.fastq >> log.txt
echo -n "FW clipped " >> report.txt && wc -l sense_FW_clipped.fastq | awk '{x=$1/4; print x}' >> report.txt

cutadapt -a "$rvprimer" -j 4 -o sense_FWRV_clipped.fastq sense_FW_clipped.fastq >> log.txt
echo -n "RV clipped " >> report.txt && wc -l sense_FWRV_clipped.fastq | awk '{x=$1/4; print x}' >> report.txt


# UNIQUE SEQUENCE COUNTING 

cat sense_FWRV_clipped.fastq | sort | uniq -c | awk '{print $1 "\t" $2}' | sort -nr -k 1 > preprocessed.tsv	