import os
from Bio import SeqIO

def get_fasta_sequences(fasta_filename):
    fasta_alignment = SeqIO.parse(fasta_filename, 'fasta')
    acc2seq = {}
    for fasta in fasta_alignment:
        acc2seq[fasta.id.replace('/', '0')] = fasta.seq
    return acc2seq

acc2seq_trimmed_all = get_fasta_sequences('data/proteome-tree/sequences-all-class-filtered-andseeds.trimmed.fa')
acc2seq_trimmed_fungi = get_fasta_sequences('data/proteome-tree/sequences-fungi-order-notFiltered.trimmed.fa')
acc2seq_seeds_trimmed = get_fasta_sequences('data/seeds-trimmed.fa')

def write_fasta(ids, output_filename):
    with open(output_filename, 'w') as outfile:
        for id in ids:
            if id in acc2seq_trimmed_all:
                outfile.write(f">{id}\n{acc2seq_trimmed_all[id]}\n")
            elif id in acc2seq_trimmed_fungi:
                outfile.write(f">{id}\n{acc2seq_trimmed_fungi[id]}\n")
            elif id in acc2seq_seeds_trimmed:
                outfile.write(f">{id}\n{acc2seq_seeds_trimmed[id]}\n")

# All (for architecture figure)
dir = "data/mrbayes/0816" 
member_dir = f"{dir}/clades/members"
clades = os.listdir(member_dir)
out_dir = f"{dir}/clades/fastas" 
for clade in clades:
    if clade.endswith('.DS_Store'):
        continue
    with open(f'{member_dir}/{clade}') as f:
        accs = f.read().splitlines()
    write_fasta(accs, f"{out_dir}/{clade}.fa")

# Fungi
dir = "data/mrbayes/short-fungal-0507" 
member_dir = f"{dir}/clades/members"
clades = os.listdir(member_dir)
out_dir = f"{dir}/clades/fastas" 
for clade in clades:
    if clade.endswith('.DS_Store'):
        continue
    with open(f'{member_dir}/{clade}') as f:
        accs = f.read().splitlines()
    write_fasta(accs, f"{out_dir}/{clade}.fa")