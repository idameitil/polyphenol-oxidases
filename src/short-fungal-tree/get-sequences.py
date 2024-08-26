from Bio import SeqIO

# read accs
acc_file = 'data/short-fungal-tree/accs'
with open(acc_file) as file:
    accs = [line.rstrip() for line in file]

# read fasta
fasta = 'data/epa-ng/fungi-order/ref-query.fa'
fasta_sequences = SeqIO.parse(fasta, 'fasta')

# write
with open('data/short-fungal-tree/short-fungal.fasta', 'w') as outfile:
    for fasta in fasta_sequences:
        if fasta.id in accs:
            outfile.write(f">{fasta.id}\n{fasta.seq}\n")