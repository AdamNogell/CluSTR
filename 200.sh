#!/bin/bash

# $1 = FWreads_file.fastq
# $2 = primers_file.txt

touch report.txt
touch log.txt

echo -n "input " >> report.txt && wc -l "$1" | awk '{x=$1/4; print x}' >> report.txt

fwprimer=$(cat "$2" | cut -f 2 | sed -n '1p')
cutadapt -g "$fwprimer" -j 4 -o FWreads_clipped.fastq "$1" >> log.txt
echo -n "clipped " >> report.txt && wc -l FWreads_clipped.fastq | awk '{x=$1/4; print x}' >> report.txt

cat FWreads_clipped.fastq | sort | uniq -c | awk '{print $1 "\t" $2}' | sort -nr -k 1 | head -n 200 > FWreads_200.tsv