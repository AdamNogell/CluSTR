#!/usr/bin/env python3

from Bio import SeqIO
import sys

# Extract centroid sequence IDs from CD-HIT-EST output file
centroid_ids = {}
with open(sys.argv[1], "r") as f:
    for line in f:
        if line.startswith(">Cluster"):
            continue  # Skip cluster header lines
        if "*" in line:
            centroid_id = line.strip().split()[2].rstrip("...").lstrip(">")
            centroid_ids[centroid_id] = centroid_id

# Extract centroid sequences from input FASTA file
output_sequences = {}
with open(sys.argv[2], "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_id = record.id.split()[0]  # Remove any additional characters after the ID
        if fasta_id in centroid_ids:
            output_sequences[fasta_id] = record.seq

# Write extracted sequences to a new FASTA file
with open(sys.argv[3], "w") as extracted:
    for seq_id, seq in output_sequences.items():
        extracted.write(f">{seq_id}\n{seq}\n")
