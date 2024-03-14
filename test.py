#!/usr/bin/env python3

from Bio import SeqIO

# Extract centroid sequences from CD-HIT-EST output file
clusters = {}
with open("cd-hit-est_output.clstr", "r") as f:
    current_cluster_id = None
    current_centroid = None
    for line in f:
        if line.startswith(">Cluster"):
            current_cluster_id = line.strip().split()[-1]
            clusters[current_cluster_id] = []
        elif "*" in line:
            clusters[current_cluster_id].append(line.strip().split()[2])

# Extract centroid sequences from input FASTA file
output_sequences = {}
with open("input.fasta", "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        for cluster_id, centroid_sequences in clusters.items():
            if record.id in centroid_sequences:
                output_sequences[record.id] = record.seq

# Write extracted sequences to a new FASTA file
with open("representative_sequences.fasta", "w") as output_file:
    for seq_id, seq in output_sequences.items():
        output_file.write(f">{seq_id}\n{seq}\n")

