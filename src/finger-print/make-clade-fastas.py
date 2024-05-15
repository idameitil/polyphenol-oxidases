import os
from Bio import SeqIO

def get_fasta_sequences(fasta_filename):
    fasta_alignment = SeqIO.parse(fasta_filename, 'fasta')
    acc2seq = {}
    for fasta in fasta_alignment:
        acc2seq[fasta.id.split('.')[0]] = fasta.seq
    return acc2seq

acc2seq_trimmed = get_fasta_sequences('data/pfam/PF00264-trimmed.fa')
acc2seq_seeds_trimmed = get_fasta_sequences('data/seeds-trimmed.fa')

def write_fasta(ids, output_filename):
    with open(output_filename, 'w') as outfile:
        for id in ids:
            if id in acc2seq_trimmed:
                outfile.write(f">{id}\n{acc2seq_trimmed[id]}\n")
            elif id in acc2seq_seeds_trimmed:
                outfile.write(f">{id}\n{acc2seq_seeds_trimmed[id]}\n")
                
dir = "data/mrbayes/all/clades/members"
clades = os.listdir(dir)
out_dir = "data/mrbayes/all/clades/fastas" 
for clade in clades:
    with open(f'{dir}/{clade}') as f:
        accs = f.read().splitlines()
    write_fasta(accs, f"{out_dir}/{clade}.fa")