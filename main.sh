#!/bin/bash

# $1 = reads1.fq.gz
# $2 = reads2.fq.gz
# $3 = primers.txt
# $4 = adapters_regions.txt
# $5 = reference_genome.fasta


#CREATING NEEDED DIRECTORIES, FILES AND VARIABLES
mkdir split mapping clustering
touch report.txt log.txt clustering/clustering_log.txt
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

echo "Unzipping completed (1/8)"

#ADAPTER TRIMMING
cutadapt -a "$adapter_3" -A "$adapter_5" -m 60 -q 18 -j 4 -o split/file_1_trimmed.fastq -p split/file_2_trimmed.fastq input1.fastq input2.fastq >> log.txt
echo -n "file 1 after adapter cutting: " >> report.txt && 
    wc -l split/file_1_trimmed.fastq | 
    awk '{print int($1/4)}' >> report.txt
echo -n "file 2 after adapter cutting: " >> report.txt && 
    wc -l split/file_2_trimmed.fastq | 
    awk '{print int($1/4)}' >> report.txt

echo "Adapter trimming completed (2/8)"

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

echo "Split completed (3/8)"

#removing unnecessary files
rm input*.fastq split/file_1_trimmed-*.fastq split/file_2_trimmed-*.fastq

#FASTQC - UNHASH IF NEEDED
#fastqc FWreads.fastq
#fastqc RVreads.fastq

#PRIMER CUTTING
cutadapt -g "$rv_primer" -j 4 -o RVreads_primerless.fastq split/RVreads.fastq >> log.txt
echo -n "RV reads after primer cut: " >> report.txt && 
    wc -l RVreads_primerless.fastq | 
    awk '{print int($1/4)}' >> report.txt

echo "Primer cutting completed (4/8)"

#MAPPING
bwa index "$5" 2> log.txt
bwa mem -t 4 "$5" RVreads_primerless.fastq > mapping/RVreads.sam 2> log.txt
samtools view -b mapping/RVreads.sam > mapping/RVreads.bam && 
    samtools sort mapping/RVreads.bam -o mapping/RVreads_sorted.bam 

echo "Mapping completed (5/8)"

#FILTERING
#filtering reads that were clipped from the left side and and pass throught the whole region of interest (FURTHER PROCESSING REQIRED)
samtools view -H mapping/RVreads_sorted.bam > header.sam && 
    samtools view mapping/RVreads_sorted.bam | 
    awk -v start="$start_pos" -v end="$end_pos" '{if ($4 <= start && $4 + length($10) >= end && $6 ~ /^[0-9]+S/) {split($6, arr, "S"); if ($4 <= start && ($4 + length($10) - arr[1]) >= end) print}}' | 
    cat header.sam - | 
    samtools view -b > filtered1.bam && 
        rm header.sam

#filtering reads that were not clipped from the left side (READY FOR MERGING)
samtools view -H mapping/RVreads_sorted.bam > header.sam && 
    samtools view mapping/RVreads_sorted.bam | 
    awk -v start="$start_pos" -v end="$end_pos" '{if ($4 <= start && $4 + length($10) >= end && $6 !~ /^[0-9]+S/) print}' | 
    cat header.sam - | 
    samtools view -b > filtered2.bam && 
        rm header.sam

#filtering reads 2 that were not clipped from the right side (READY FOR MERGING)
samtools view -H filtered2.bam > header.sam && 
    samtools view filtered2.bam | 
    awk '{if ($6 !~ /S$/) print}' | 
    cat header.sam - | 
    samtools view -b > filtered2A.bam && 
        rm header.sam

#filtering reads 2 that were clipped from the right side (FURTHER PROCESSING REQUIRED)
samtools view -H filtered2.bam > header.sam && 
    samtools view filtered2.bam | 
    awk '{if ($6 ~ /S$/) print}' | 
    cat header.sam - | 
    samtools view -b > filtered2B.bam && 
        rm header.sam

#filtering reads 2B that passed through the whole region of interest (READY FOR MERGING)
samtools view -H filtered2B.bam > header.sam && 
    samtools view filtered2B.bam | 
    awk -v start="$start_pos" -v end="$end_pos" '{split($6, arr, "S"); n=split(arr[1], ar, "M"); if ($4 <= start && ($4 + $10 - ar[n] >= end)) print}' | 
    cat header.sam - | 
    samtools view -b > filtered2C.bam && 
        rm header.sam

#filtering reads that were not clipped from the right side (READY FOR MERGING)
samtools view -H filtered1.bam > header1.sam && 
    samtools view filtered1.bam | 
    awk '{if ($6 !~ /S$/) print}' | 
    cat header1.sam - | 
    samtools view -b > filtered3.bam && 
        rm header1.sam

#filtering reads that were clipped from the right side (FURTHER PROCESSING REQUIRED)
samtools view -H filtered1.bam > header1.sam && 
    samtools view filtered1.bam | 
    awk '{if ($6 ~ /S$/) print}' | 
    cat header1.sam - | 
    samtools view -b > filtered4.bam && 
        rm header1.sam

#filtering reads that contains deletion (FURTHER PROCESSING REQUIRED)
samtools view -H filtered4.bam > header1.sam && 
    samtools view filtered4.bam | 
    awk '{if ($6 ~ /D/) print}' | 
    cat header1.sam - | 
    samtools view -b > filtered5.bam && 
        rm header1.sam

#filtering reads that do not contain deletion (READY FOR MERGING)
samtools view -H filtered4.bam > header1.sam && 
    samtools view filtered4.bam | 
    awk '{if ($6 !~ /D/) print}' | 
    cat header1.sam - | 
    samtools view -b > filtered6.bam && 
        rm header1.sam

#filtering reads 5 that passed through the whole region of interest (READY FOR MERGING)
samtools view -H filtered5.bam > header1.sam && 
    samtools view filtered5.bam | 
    awk -v start="$start_pos" -v end="$end_pos" '{split($6, arr, "S"); split(arr[2], ar, "M"); if ($4 <= start && ($4 + length($10) - arr[1] - ar[3]) >= end) print}' | 
    cat header1.sam - | 
    samtools view -b > filtered7.bam && 
        rm header1.sam

#filtering reads 6 that passed through the whole region of interest (READY FOR MERGING)
samtools view -H filtered6.bam > header1.sam && 
    samtools view filtered6.bam | 
    awk -v start="$start_pos" -v end="$end_pos" '{split($6, arr, "S"); split(arr[2], ar, "M"); if ($4 <= start && ($4 + length($10) - arr[1] - ar[2]) >= end) print}' | 
    cat header1.sam - | 
    samtools view -b > filtered8.bam && 
        rm header1.sam

#merging the correct files into a single one
samtools merge -h  mapping/RVreads_sorted.bam -o mapping/preclustering_mapped.bam filtered2A.bam filtered2C.bam filtered3.bam filtered7.bam filtered8.bam && 
    samtools index mapping/preclustering_mapped.bam

#converting bam to fasta
samtools fasta mapping/preclustering_mapped.bam > filtered.fasta

echo "Filtering completed (6/8)"

#CLUSTERING
cd-hit-est -i filtered.fasta -o clustering/cluster_c96_n10 -T 0 -M 0 -n 10 -c 0.96 -d 0 >> clustering/clustering_log.txt
mv clustering/cluster_c96_n10 clustering/cluster_c96_n10.fasta

echo "Clustering completed (7/8)"

#extracting representative sequences (centroids) from each cluster, then cutting the reads to desired region, then filtering frequently occuring sequences
../scripts/main.py "$4"

echo "Python script successfully executed - said Bash (8/8)"

#removing unnecessary files
rm filtered1.bam filtered2.bam filtered2A.bam filtered2B.bam filtered2C.bam filtered3.bam filtered4.bam filtered5.bam filtered6.bam filtered7.bam filtered8.bam
rm cigar.fasta cluster_count.txt cut.fasta extracted.fasta mapped* final_temp* filtered.fasta RVreads_primerless.fastq
rm mapping/RV*

echo "hotovo:)"