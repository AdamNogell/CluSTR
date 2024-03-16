#!/usr/bin/env python3

# argv[1] = clustered .clstr file
# argv[2] = clustered .fasta file
# argv[3] = output .txt file

from Bio import SeqIO
import sys

# Extract centroid sequence IDs from CD-HIT-EST output file
centroid_ids = {}
with open(sys.argv[1], "r") as f:
    for line in f:
        if line.startswith(">Cluster"):
            current_cluster_id = line.strip().split()[-1]
            continue
        if "*" in line:
            centroid_id = line.strip().split()[2].rstrip("...").lstrip(">")
            centroid_ids[centroid_id] = current_cluster_id

# Extract centroid sequences from input FASTA file
output_sequences = {}
with open(sys.argv[2], "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_id = record.id.split()[0]
        if fasta_id in centroid_ids:
            output_sequences[centroid_ids[fasta_id]] = record.seq

# Write extracted sequences to a new FASTA file
with open(sys.argv[3], "w") as extracted:
    sorted_sequences = sorted(output_sequences.items(), key=lambda x: int(x[0]))
    for cluster_id, seq in sorted_sequences:
        extracted.write(f">{cluster_id}\n{seq}\n")
