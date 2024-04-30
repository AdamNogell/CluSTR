#!/usr/bin/env python3

from Bio import SeqIO
import re, sys, subprocess

def parseCigar(C:str): 
    """
    Parse CIGAR string into a list of tuples
    
    C:           GIGAR
    """
    Y = re.findall(r'(\d+)([MIDNSHP=X])', C)
    Y = [(int(y[0]), y[1]) for y in Y]
    return Y

def applyCigar(start_pos:int, cgr:str, sq:str):
    """
    Rewrite the sequence according to the CIGAR and the starting position
    
    start_pos:     Starting position of the alignment
    cgr:           GIGAR
    sq:            Nucleotide sequence
    """
    digits = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    aln_sq = ''
    current_position = 0
    operation_count = ''
    parsed = parseCigar(cgr)
    found_S = False
    for item in parsed:
        if 'S' in item:
            found_S = True
            S_item = item
            for x in S_item:
                if isinstance(x, int):
                    S = x
            break

    if found_S:
        aln_sq += '-' * (int(start_pos) - 1 - (S))
    else:
        aln_sq += '-' * (int(start_pos) - 1)

    for c in cgr:
        if c in digits:
            operation_count += c
        else:
            operation_count = int(operation_count)
            if c == 'M' or c == 'I':
                aln_sq += sq[current_position : current_position + operation_count]
                current_position += operation_count
            elif c == 'D':
                aln_sq += '-' * operation_count
            elif c == 'S' or c == 'H' or c == 'N':
                aln_sq += '*' * operation_count
                current_position += operation_count
            else:
                return None
            
            operation_count = ''

    return aln_sq

with open(sys.argv[1], 'r') as region:
    line = region.readlines()
    lines = [x.strip() for x in line]
    start_pos = int(lines[2])
    end_pos = int(lines[3])

for record in SeqIO.parse('../U13369/U13369.1.fasta', 'fasta'):
    ref_seq = str(record.seq[start_pos:end_pos])

# Extract centroid sequence IDs from CD-HIT-EST output file
centroid_ids = {}
with open('clustering/cluster_c96_n10.clstr', "r") as f:
    for line in f:
        if line.startswith(">Cluster"):
            current_cluster_id = line.strip().split()[-1]
            continue
        if "*" in line:
            centroid_id = line.strip().split()[2].rstrip("...").lstrip(">")
            centroid_ids[centroid_id] = current_cluster_id

# Extract centroid sequences from input FASTA file
output_sequences = {}
with open('clustering/cluster_c96_n10.fasta', "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_id = record.id.split()[0]
        if fasta_id in centroid_ids:
            output_sequences[centroid_ids[fasta_id]] = record.seq

# Write extracted sequences to a new FASTA file
with open('extracted.fasta', "w") as extracted:
    sorted_sequences = sorted(output_sequences.items(), key=lambda x: int(x[0]))
    for cluster_id, seq in sorted_sequences:
        extracted.write(f">{cluster_id}\n{seq}\n")

subprocess.run('bwa mem -t 4 ../U13369/U13369.1.fasta extracted.fasta > mapped.sam 2> clustering/clustering_log.txt', shell=True)
subprocess.run('samtools view -b mapped.sam > mapped.bam && samtools sort mapped.bam > mapping/postclustering_mapped.bam && samtools index mapping/postclustering_mapped.bam', shell=True)

"""
Input is a .sam file that is parsed according to columns. Then columns 4 (starting position), 6 (CIGAR) and 10 (sequence) are passed down to the applyCigar function that outputs a sequence cut according to the coordinates of the region of interest. 
Deletions are replaced in the output with '-' and sequences that have been hard-clipped or soft-clipped are replaced with '*'. The output is .fasta file with the corresponding sequence ID and then the modified sequence in the next line.
"""
split = []
with open('mapped.sam', "r") as sam, open('cigar.fasta', "w") as fasta:
    for line in sam:
        if '@' not in line:
            split = line.split()
            parsed = applyCigar(split[3], split[5], split[9])
            fasta.write(f'>{split[0]}\n{parsed[start_pos:end_pos]}\n')

# Find the sequences that doesn't contain '*' in the last .fasta output file and rewrite them into new .fasta file with corresponding sequence ID
with open('cigar.fasta', 'r') as cigar, open('cut.fasta', 'w') as cut:
    lines = cigar.readlines()
    for i in range(1, len(lines)):
        if not any(c in lines[i] for c in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '*']):
            cut.write(lines[i-1])
            cut.write(lines[i])

# Read the FASTA file and store the sequences in a dictionary
sequences = {}
for record in SeqIO.parse('cut.fasta', 'fasta'):
    sequence_id = record.id
    sequence = str(record.seq)
    if sequence not in sequences:
        sequences[sequence] = set()
    sequences[sequence].add(sequence_id)

#Modify the dictionary into [[list of [lists]] of [string and [list]]]
sequence_values =[[sequence, list(sorted(value))] for sequence, value in sequences.items()] # [['seq', ['cluster#']], ['seq', ['cluster#']]]

#Count sequences for each cluster and output them into cluster_count.txt file
cluster_count = 0
with open('clustering/cluster_c96_n10.clstr', "r") as A, open('cluster_count.txt', "w") as B:
    for line in A:
        if line.startswith(">Cluster"):
            line = line.split()
            if cluster_count > 0:
                B.write(f"{cluster_count}\n")
                B.write(f"{line[1]}\t")
                cluster_count = 0
            else:
                B.write(f"{line[1]}\t")
                cluster_count = 0
        else:
            cluster_count += 1
    if cluster_count > 0:
        B.write(f"{cluster_count}")

#Combine cluster_count and sequence_values into a single output that tells how many unique sequences occur in the original file
with open('cluster_count.txt', 'r') as A, open('final_temp.txt', 'w') as B:
    Alines = A.readlines()
    for cluster in sequence_values:
        temp = 0
        for x in cluster[1]: 
            for line in Alines:
                linespl = line.split()
                if x in linespl[0]:
                    temp += int(linespl[1])
        cluster[1] = temp
    for final in sequence_values:
        a = str(final[0])
        b = str(final[1])
        B.write(f"{a}\t{b}\n")

subprocess.run('sort -k2nr final_temp.txt | awk \'{if ($2 >= 1000) print}\' > final_temp2.txt', shell=True)

with open('final_temp2.txt', 'r') as temp, open('final.txt', 'w') as final:
    final.write(f"{ref_seq}\tREFERENCE\n")
    for line in temp:
        final.write(line)


print("Python script successfully executed - said Python (8/8)")