#!/bin/bash

# $1 = reads1.fq.gz
# $2 = reads2.fq.gz
# $3 = primers.txt
# $4 = adapters_regions.txt
# $5 = reference_genome.fasta


#CREATING NEEDED DIRECTORIES AND FILES
mkdir split
mkdir top
mkdir mapping
mkdir clustering
touch report.txt
touch log.txt
touch clustering/log_c96_n10.txt
touch clustering/extracted_c96_n10.fasta
rv_primer=$( cut -f 2 "$3" | sed -n '2p')
adapter_5=$(sed -n '1p' "$4")
adapter_3=$(sed -n '2p' "$4")
start_pos=$(sed -n '3p' "$4")
end_pos=$(sed -n '4p' "$4")
echo "Number of reads" >> report.txt

#UNZIPPING READS
zcat -f "$1" > input1.fastq
echo -n "input file 1: " >> report.txt &&
    wc -l input1.fastq | 
    awk '{print int($1/4)}' >> report.txt

zcat -f "$2" > input2.fastq
echo -n "input file 2: " >> report.txt &&
    wc -l input2.fastq | 
    awk '{print int($1/4)}' >> report.txt

#ADAPTER TRIMMING
cutadapt -a "$adapter_3" -A "$adapter_5" -m 60 -q 18 -j 4 -o split/file_1_trimmed.fastq -p split/file_2_trimmed.fastq input1.fastq input2.fastq >> log.txt
echo -n "file 1 after adapter cutting: " >> report.txt && 
    wc -l split/file_1_trimmed.fastq | 
    awk '{print int($1/4)}' >> report.txt
echo -n "file 2 after adapter cutting: " >> report.txt && 
    wc -l split/file_2_trimmed.fastq | 
    awk '{print int($1/4)}' >> report.txt

#SPLITTING FW AND RV READS
../fastx_barcode_splitter.pl --bcfile "$3" --bol --prefix split/file_1_trimmed- --suffix .fastq --mismatches 5 --partial 5 < split/file_1_trimmed.fastq >> log.txt
echo -n "file 1 FW reads: " >> report.txt && 
    wc -l split/file_1_trimmed-FW.fastq | 
    awk '{print int($1/4)}' >> report.txt
echo -n "file 1 RV reads: " >> report.txt && 
    wc -l split/file_1_trimmed-RV.fastq | 
    awk '{print int($1/4)}' >> report.txt
../fastx_barcode_splitter.pl --bcfile "$3" --bol --prefix split/file_2_trimmed- --suffix .fastq --mismatches 5 --partial 5 < split/file_2_trimmed.fastq >> log.txt
echo -n "file 2 FW reads: " >> report.txt && 
    wc -l split/file_2_trimmed-FW.fastq | 
    awk '{print int($1/4)}' >> report.txt
echo -n "file 2 RV reads: " >> report.txt && 
    wc -l split/file_2_trimmed-RV.fastq | 
    awk '{print int($1/4)}' >> report.txt

#joining FW files together and RV files together
cat split/file_1_trimmed-FW.fastq split/file_2_trimmed-FW.fastq > split/FWreads.fastq
echo -n "FW reads total: " >> report.txt && 
    wc -l split/FWreads.fastq | 
    awk '{print int($1/4)}' >> report.txt
cat split/file_1_trimmed-RV.fastq split/file_2_trimmed-RV.fastq > split/RVreads.fastq
echo -n "RV reads total: " >> report.txt && 
    wc -l split/RVreads.fastq | 
    awk '{print int($1/4)}' >> report.txt

#removing unnecessary files
rm input1.fastq
rm input2.fastq
rm split/file_1_trimmed-FW.fastq
rm split/file_1_trimmed-RV.fastq
rm split/file_2_trimmed-FW.fastq
rm split/file_2_trimmed-RV.fastq

#FASTQC - UNHASH IF NEEDED
#fastqc FWreads.fastq
#fastqc RVreads.fastq

#PRIMER CUTTING
cutadapt -g "$rv_primer" -j 4 -o top/RVreads_primerless.fastq split/RVreads.fastq >> log.txt
echo -n "RV reads after primer cut: " >> report.txt && 
    wc -l top/RVreads_primerless.fastq | 
    awk '{print int($1/4)}' >> report.txt

#creating a list of 200 most frequent RV reads
awk '{NR%4 == 2}' top/RVreads_primerless.fastq > top/RVreads_primerless_awk.txt
sort top/RVreads_primerless_awk.txt | 
uniq -c | 
awk '{print $1 "\t" $2}' | 
sort -nr -k 1 | 
head -n 200 > top/RVreads_top200.tsv

#MAPPING
bwa index "$5" >> log.txt
bwa mem -t 4 "$5" top/RVreads_primerless.fastq > mapping/RVreads.sam
samtools view -b mapping/RVreads.sam > mapping/RVreads.bam && 
    samtools sort mapping/RVreads.bam -o mapping/RVreads_sorted.bam && 
        samtools index mapping/RVreads_sorted.bam
samtools view -b -o mapping/RV_region_sorted.bam mapping/RVreads_sorted.bam U13369.1:41842-41877 && 
    samtools index mapping/RV_region_sorted.bam

#FILTERING
#filtering all reads that were clipped from the left side and pass through the whole Q2
samtools view -H mapping/RVreads_sorted.bam > header1.sam && 
    samtools view mapping/RVreads_sorted.bam | 
    awk -v start="$start_pos" -v end="$end_pos" '{if ($4 <= start && $4 + length($10) >= end && $6 ~ /^[0-9]+S/) {split($6, arr, "S"); if ($4 <= start && ($4 + length($10) - arr[1]) >= end) print}}' | 
    cat header1.sam - | 
    samtools view -b > filtered1.bam && 
        rm header1.sam

#filtering all reads that were not clipped from the left side and pass through the whole Q2
samtools view -H mapping/RVreads_sorted.bam > header2.sam && 
    samtools view mapping/RVreads_sorted.bam | 
    awk -v start="$start_pos" -v end="$end_pos" '{if ($4 <= start && $4 + length($10) >= end && $6 !~ /^[0-9]+S/) print}' | 
    cat header2.sam - | 
    samtools view -b > filtered2.bam && 
        rm header2.sam

#merging all reads that pass through Q2 and indexe the merged bam file
samtools merge mapping/filtered.bam filtered1.bam filtered2.bam && 
    samtools index mapping/filtered.bam

#converting bam to fasta
samtools fasta mapping/filtered.bam > clustering/filtered.fasta

#removing unnecessary files
rm filtered1.bam
rm filtered2.bam

#CLUSTERING
cd-hit-est -i clustering/filtered.fasta -o clustering/cluster_c96_n10 -T 0 -M 0 -n 10 -c 0.96 -d 0 >> clustering/log_c96_n10.txt
mv clustering/cluster_c96_n10 clustering/cluster_c96_n10.fasta

#extracting representative sequences (centroids) from each cluster
../scripts/cluster_extract.py clustering/cluster_c96_n10.clstr clustering/cluster_c96_n10.fasta clustering/extracted_c96_n10.fasta

#mapping the centroids
bwa mem -t 4 "$5" clustering/extracted_c96_n10.fasta > clustering/aln_c96_n10.sam
samtools view -b clustering/aln_c96_n10.sam > clustering/aln_c96_n10.bam &&
    samtools sort clustering/aln_c96_n10.bam > clustering/aln_c96_n10_sorted.bam &&
        samtools index clustering/aln_c96_n10_sorted.bam
