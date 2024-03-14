#!/bin/bash

# Extract the first 30 bases from each read
awk 'NR%4 == 2 {print substr($0, 1, 10)}' "$1" > extracted_sequences.txt

# Count occurrences of each sequence, sort by count in descending order, and save to output file
cat extracted_sequences.txt | sort | uniq -c | awk '{print $1 "\t" $2}' | sort -nr -k 1 | head -n 200 > output.tsv
