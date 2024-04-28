#!/bin/bash

# $1 = reads1.fq.gz
# $2 = reads2.fq.gz
# $3 = primers.txt
# $4 = adapters_regions.txt
# $5 = reference_genome.fasta


#CREATING NEEDED DIRECTORIES, FILES AND VARIABLES
mkdir split top mapping clustering
touch report.txt log.txt clustering/log_c96_n10.txt clustering/extracted_c96_n10.fasta region.bed
rv_primer=$( cut -f 2 "$3" | sed -n '2p')
adapter_5=$(sed -n '1p' "$4")
adapter_3=$(sed -n '2p' "$4")
start_pos=$(sed -n '3p' "$4")
declare -i x="$start_pos"
region_pos=$((x - 1))
end_pos=$(sed -n '4p' "$4")
printf "U13369.1\t0\t%d\nU13369.1\t%d\t43999" "$region_pos" "$end_pos" > region.bed
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
rm input1.fastq input2.fastq split/file_1_trimmed-FW.fastq split/file_1_trimmed-RV.fastq split/file_2_trimmed-FW.fastq split/file_2_trimmed-RV.fastq

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
samtools merge -h  mapping/RVreads_sorted.bam -o mapping/filtered.bam filtered2A.bam filtered2C.bam filtered3.bam filtered7.bam filtered8.bam && 
    samtools index mapping/filtered.bam

#removing unnecessary files
rm filtered1.bam filtered2.bam filtered2A.bam filtered2B.bam filtered2C.bam filtered3.bam filtered4.bam filtered5.bam filtered6.bam filtered7.bam filtered8.bam

#REGION CUT
#cutting all reads to contain only the region of interest
samtools ampliconclip -o clustering/region_cut.bam --hard-clip --both-ends -b region.bed mapping/filtered.bam &&
    samtools sort clustering/region_cut.bam > clustering/region_cut_sorted.bam &&
rm region.bed

#converting bam to fasta
samtools fasta clustering/region_cut_sorted.bam > clustering/region_cut_sorted.fasta

#CLUSTERING
cd-hit-est -i clustering/region_cut_sorted.fasta -o clustering/cluster_c96_n10 -T 0 -M 0 -n 10 -c 0.96 -d 0 >> clustering/log_c96_n10.txt
mv clustering/cluster_c96_n10 clustering/cluster_c96_n10.fasta

#extracting representative sequences (centroids) from each cluster
../scripts/cluster_extract.py clustering/cluster_c96_n10.clstr clustering/cluster_c96_n10.fasta clustering/extracted_c96_n10.fasta
