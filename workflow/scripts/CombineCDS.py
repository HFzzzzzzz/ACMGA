#!python3
from Bio import SeqIO
import os
seen_sequences = set()
directory_path = './CDS'

all_files = os.listdir(directory_path)

fasta_files = [file for file in all_files if file.endswith("fasta")]
# print(fasta_files)
with open("non_duplicate_CDS.fa", "w") as output_handle:
    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.seq not in seen_sequences:
                seen_sequences.add(record.seq)
                print(">"+record.name)
                print(record.seq)
