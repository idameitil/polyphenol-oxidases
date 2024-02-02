from Bio import SeqIO

# Read filtered
filtered_fasta_filename = 'data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta'
fasta_sequences = SeqIO.parse(filtered_fasta_filename, 'fasta')
accessions_include = []
for fasta in fasta_sequences:
    accessions_include.append(fasta.id)

# Read alignment
alignment_filename = 'data/pfam/PF00264.alignment.uniprot-cleaned.fa'
fasta_sequences = SeqIO.parse(alignment_filename, 'fasta')

# Write output
output_filename = 'data/pfam/PF00264.alignment.uniprot-cleaned-filtered.fa'
accessions_done = []
with open(output_filename, 'w') as outfile:
    for fasta in fasta_sequences:
        # print(fasta.id.split('.')[0])
        if fasta.id.split('.')[0] in accessions_include and fasta.id.split('.')[0] not in accessions_done:
            outfile.write(f">{fasta.id.split('.')[0]}\n{fasta.seq}\n")
            accessions_done.append(fasta.id.split('.')[0])