import pandas as pd
import string

table = str.maketrans('', '', string.ascii_lowercase)

data = {}
input_filename = 'data/pfam/PF00264.alignment.uniprot'
output_filename = 'data/pfam/PF00264.alignment.uniprot-cleaned.fa'
count = 0
with open(input_filename, 'r') as infile:
    with open(output_filename, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue
            else:
                id = line[:31].strip()
                seq = line[31:].strip()
                translated_seq = seq.translate(table)
                new_seq = translated_seq.replace('.', '')
                outfile.write(f">{id}\n{new_seq}\n")
