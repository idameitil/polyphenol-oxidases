from Bio import SeqIO

# read accs
acc_file = 'data/short-fungal-tree/accs'
with open(acc_file) as file:
    accs = [line.rstrip() for line in file]

# read fasta
fasta = 'data/epa-ng/fungi-order/ref-query.fa'
fasta_sequences = SeqIO.parse(fasta, 'fasta')

# read seed fasta
seed_fasta = 'data/seeds-trimmed.fa'
seed_sequences = SeqIO.parse(seed_fasta, 'fasta')

seeds_include = ['AoCO4', 'CgPPO-473', 'CgPPO-266', 'MtPPO-809', 'MtPPO-010', 'PpPPO-c2092', 'MtPPO7']

# write
with open('data/short-fungal-tree/short-fungal.fasta', 'w') as outfile:
    for fasta in fasta_sequences:
        if fasta.id in accs:
            outfile.write(f">{fasta.id}\n{fasta.seq}\n")
    for fasta in seed_sequences:
        if fasta.id in seeds_include:
            outfile.write(f">{fasta.id}\n{fasta.seq}\n")