from Bio import SeqIO

fasta_sequences = SeqIO.parse('data/pfam/protein-matching-PF00264-shortheaders.fasta', 'fasta')

output_filename = 'data/pfam/protein-matching-PF00264-shortheaders-70-1000.fasta'

with open(output_filename, 'w') as outfile:
    for fasta in fasta_sequences:
        if len(fasta.seq) > 70 and len(fasta.seq) < 1000:
            outfile.write(f">{fasta.id}\n{fasta.seq}\n")